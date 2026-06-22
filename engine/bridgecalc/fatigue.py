"""疲勞驗核（P1）：鋼腱應力幅、混凝土壓疲勞、箍筋疲勞。

對應公式卡 P1（疲勞設計）、算例_疲勞設計驗核、疲勞載重計算快查卡。
疲勞載重 = 疲勞卡車（35/145/145、後軸固定 9m、單車道）× IM 15%；Fatigue I γ=1.75。
單位：力 N、長度 mm、應力 MPa；彎矩/剪力介面 kN·m / kN。
"""
from dataclasses import dataclass
from .model import Section


@dataclass
class FatigueResult:
    dsig_ps: float    # 鋼腱應力幅 MPa
    sig_c_max: float  # 混凝土最大壓應力（疲勞）MPa
    ps_ok: bool       # Δσ_ps ≤ 125？
    c_ok: bool        # σ_c ≤ 0.40 f'c？
    CR_ps: float
    CR_c: float


def fatigue_check(section: Section, Pe: float, e: float,
                  M_perm_kNm: float, dM_fatigue_kNm: float, fc: float,
                  Ep_Ec: float = 6.6, gamma: float = 1.75) -> FatigueResult:
    """鋼腱應力幅與混凝土壓疲勞。

    M_perm：永久載重彎矩（DC+SDL）；dM_fatigue：疲勞卡車彎矩幅（含 IM）。
    """
    dM = gamma * dM_fatigue_kNm * 1e6        # N·mm（含 Fatigue I 因子）
    dsig_ps = Ep_Ec * (dM / section.I) * e   # 鋼腱（底）應力幅 MPa
    # 永久載重頂緣壓應力（壓為負）
    sig_perm_top = -Pe / section.A + Pe * e / section.St - M_perm_kNm * 1e6 / section.St
    dsig_c_top = dM / section.I * section.yt  # 疲勞彎矩在頂緣增量（壓）
    sig_c_max = abs(sig_perm_top) + dsig_c_top
    limit_c = 0.40 * fc
    return FatigueResult(dsig_ps, sig_c_max,
                         dsig_ps <= 125.0, sig_c_max <= limit_c,
                         dsig_ps / 125.0, sig_c_max / limit_c)


def stirrup_fatigue(dV_kN: float, s: float, Av: float, dv: float,
                    limit: float = 165.0):
    """箍筋疲勞應力幅 Δf_sv = ΔV·s/(Av·dv)（MPa）。

    dV：疲勞卡車剪力幅（kN）；s 間距、Av 箍筋面積、dv 有效剪力深（mm）。
    回傳 (Δf_sv, ok)。近支承段常控制 → 加密箍筋。
    """
    dfsv = dV_kN * 1e3 * s / (Av * dv)
    return dfsv, dfsv <= limit
