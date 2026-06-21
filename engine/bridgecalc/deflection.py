"""撓度與長期下撓／預拱（C2/C3）。

對應公式卡 C2、算例_預拱設計（C3）。
單位：L mm、Ec MPa、I mm⁴、w kN/m（= N/mm）、Pe N、e mm；撓度 mm（向下為正）。
關鍵：潛變同時作用於自重與預力（皆 ×(1+φ)）。預力反拱 w_eq = 8·Pe·a/L²。
"""
from dataclasses import dataclass
from .model import Section


@dataclass
class DeflectionResult:
    K: float          # 撓度係數 mm/(N/mm)
    d_DL: float       # 自重+SDL 彈性下撓 mm
    w_eq: float       # 預力等效上撐 kN/m
    d_PT: float       # 預力彈性上拱 mm
    net_elastic: float    # 淨彈性（>0 下撓）mm
    LBR: float        # 荷重平衡率 w_eq/w_DL
    net_long_term: float  # 淨長期（×(1+φ)）mm
    d_LL: float       # 活載即時撓度 mm
    camber: float     # 建議預拱 mm（淨長期下撓 + 沉陷）
    d_LL_ok: bool     # δ_LL ≤ L/800？


def deflection_analysis(L: float, Ec: float, section: Section,
                        w_DL: float, Pe: float, e: float, w_LL: float,
                        phi: float = 2.0, settlement: float = 5.0) -> DeflectionResult:
    K = 5 * L**4 / (384 * Ec * section.I)         # mm/(N/mm)
    d_DL = K * w_DL
    w_eq = 8 * Pe * e / L**2                        # N/mm（= kN/m）
    d_PT = K * w_eq
    net_el = d_DL - d_PT
    LBR = d_PT / d_DL
    net_LT = net_el * (1 + phi)
    d_LL = K * w_LL
    camber = max(net_LT, 0.0) + settlement
    return DeflectionResult(K, d_DL, w_eq, d_PT, net_el, LBR, net_LT, d_LL,
                            camber, d_LL <= L / 800)
