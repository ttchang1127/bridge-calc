"""腹板抗剪（D1）：軸壓 fpc、鋼腱垂直分量 Vp、莫耳圓主拉應力、Vcw、所需箍筋。

對應公式卡 D1（腹板剪力設計）、算例_腹板抗剪設計。單位：力 N、長度 mm、應力 MPa。
台灣主拉/Vcw 係數 0.3√f'c(kgf/cm²) 已換算為 MPa：0.094·√f'c（f'c MPa）。
"""
from dataclasses import dataclass
from math import sqrt
from .model import Section


@dataclass
class ShearResult:
    fpc: float        # 軸向壓應力 MPa（正值）
    slope: float      # 控制斷面鋼腱斜率（≈sinθ）
    Vp: float         # 鋼腱垂直分量（每腹板）N
    tau: float        # 剪應力 MPa
    sigma1: float     # 主拉應力 MPa
    sigma1_limit: float   # 台灣主拉容許 MPa
    sigma1_ok: bool       # σ1 ≤ 限值？（近支承常超→靠箍筋）
    Vcw: float        # 混凝土抗剪（含 Vp）N
    Vs_req: float     # 所需鋼材抗剪 N
    Av_s_req: float   # 所需箍筋 mm²/mm


def principal_tension_limit_TW(fc: float) -> float:
    """台灣主拉應力容許 = 0.3√f'c(kgf/cm²)，換算 MPa = 0.094·√f'c。"""
    return 0.094 * sqrt(fc)


def shear_web(Pe: float, section: Section, e: float, fc: float,
              bw_eff: float, dv: float, Vu: float,
              x_control: float, L: float,
              n_webs: int = 2, phi: float = 0.85, fsy: float = 420.0) -> ShearResult:
    """近支承控制斷面腹板抗剪（每腹板）。

    Pe 全斷面 N、e 跨中偏心 mm、Vu 每腹板設計剪力 N、x_control 控制斷面距支承 mm。
    """
    fpc = Pe / section.A
    slope = (4 * e / L) * (1 - 2 * x_control / L)        # 拋物線斜率 ≈ sinθ
    Vp = (Pe / n_webs) * slope
    tau = Vu / (bw_eff * dv)
    sigma1 = -fpc / 2 + sqrt((fpc / 2) ** 2 + tau ** 2)
    lim = principal_tension_limit_TW(fc)

    vc = 0.094 * sqrt(fc)                                 # 台灣混凝土拉貢獻 MPa
    Vcw = (vc + 0.3 * fpc) * bw_eff * dv + Vp
    Vs_req = (Vu - phi * Vcw) / phi
    Av_s_req = max(Vs_req, 0.0) / (fsy * dv)
    return ShearResult(fpc, slope, Vp, tau, sigma1, lim, sigma1 <= lim,
                       Vcw, Vs_req, Av_s_req)


def phiVn(Vcw: float, Av_s_provided: float, dv: float,
          phi: float = 0.85, fsy: float = 420.0) -> float:
    """實配箍筋下的設計抗剪 φ(Vc+Vs)。"""
    Vs = Av_s_provided * fsy * dv
    return phi * (Vcw + Vs)


def Av_s_min_TW(fc: float, bw_eff: float, fsy: float = 420.0) -> float:
    """台灣最小箍筋 max(0.2√f'c·b/fsy, 0.35·b/fsy)（mm²/mm）。"""
    return max(0.2 * sqrt(fc) * bw_eff / fsy, 0.35 * bw_eff / fsy)
