"""連續梁中墩（次彎矩 M2、T 斷面極限彎曲、墩斷面服務性）。

對應 算例_連續梁次彎矩與服務性應力驗算、算例_極限強度彎曲設計 §十二（中墩）。
40+40 兩跨連續後張箱梁。★ 中墩為控制斷面，兩項皆「失敗」：
  - T 斷面 M1：負彎矩使 NA 進腹板（c>hf）→ Mn 縮水 → CR≈0.42 嚴重不足
  - B 墩服務性：底緣壓應力 ≈20 MPa > 0.45f'c=18（1.11 倍）
單位：力 N、長度 mm、應力 MPa；M 以 kN·m 表示。
"""
from dataclasses import dataclass
from math import sqrt
from .model import Section
from .flexure import beta1


def secondary_moment(M_total_kNm: float, M1_kNm: float) -> float:
    """次彎矩 M2 = M_total − M1（超靜定連續梁；M1=ΣPe·e 一次彎矩）。"""
    return M_total_kNm - M1_kNm


def primary_moment(layers) -> float:
    """一次預力彎矩 M1 = Σ(Pe_layer · e_layer)。layers: [(Pe_kN, e_m), ...]。"""
    return sum(Pe * e for Pe, e in layers)


@dataclass
class TFlexureResult:
    c: float          # 中性軸 mm
    flanged: bool     # True=T 斷面（NA 進腹板）
    a: float          # 等效矩形深 mm
    fps: float        # 鋼腱應力 MPa
    Mn: float         # 公稱彎矩 kN·m
    eps_t: float      # 極限拉應變
    phi: float        # 強度折減
    phiMn: float      # 設計彎矩 kN·m
    CR: float         # 強度比 φMn/Mu
    ok: bool


def flexural_strength_T(Aps: float, fpu: float, fc: float,
                        b: float, hf: float, bw: float, dp: float,
                        Mu_kNm: float, dt: float = None) -> TFlexureResult:
    """T／矩形斷面極限彎曲：先試矩形，c>hf 改 T 斷面（壓力區進腹板）。

    b 受壓翼緣有效寬、hf 翼緣厚、bw 腹板寬、dp 鋼腱有效深。
    中墩負彎矩 → 受壓區在底板（窄）→ 常進腹板需 T 斷面公式。
    """
    fpy = 0.90 * fpu
    k = 2 * (1.04 - fpy / fpu)
    b1 = beta1(fc)
    c_rect = Aps * fpu / (0.85 * fc * b1 * b + k * Aps * fpu / dp)
    if c_rect <= hf:
        flanged = False
        c = c_rect
        a = b1 * c
        fps = fpu * (1 - k * c / dp)
        Mn = Aps * fps * (dp - a / 2) / 1e6
    else:
        flanged = True
        c = (Aps * fpu - 0.85 * fc * (b - bw) * hf) / (0.85 * fc * b1 * bw + k * Aps * fpu / dp)
        a = b1 * c
        fps = fpu * (1 - k * c / dp)
        Mn = (Aps * fps * (dp - a / 2)
              + 0.85 * fc * (b - bw) * hf * (a / 2 - hf / 2)) / 1e6
    dt = dt if dt is not None else dp
    eps_t = (dt - c) / c * 0.003
    if eps_t >= 0.005:
        phi = 1.0
    elif eps_t <= 0.002:
        phi = 0.75
    else:
        phi = 0.75 + 0.25 * (eps_t - 0.002) / 0.003
    phiMn = phi * Mn
    return TFlexureResult(c, flanged, a, fps, Mn, eps_t, phi, phiMn,
                          phiMn / Mu_kNm, phiMn >= Mu_kNm)


def pier_service_stress(Pe: float, section: Section, e: float, M_ext_kNm: float):
    """連續梁墩斷面服務性應力（壓為負）。e 形心上為正（負彎矩區頂板 PT）。

    σ_t = −Pe/A + Pe·e/St − M_ext/St ；σ_b = −Pe/A − Pe·e/Sb + M_ext/Sb。
    M_ext 含 M2（= M_DL+M_LL+M2）。回傳 (σ_t, σ_b)。
    """
    M = M_ext_kNm * 1e6
    st = -Pe / section.A + Pe * e / section.St - M / section.St
    sb = -Pe / section.A - Pe * e / section.Sb + M / section.Sb
    return st, sb
