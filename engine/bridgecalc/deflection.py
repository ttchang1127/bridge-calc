"""撓度與長期下撓／預拱（C2/C3）。

對應公式卡 C2、算例_預拱設計（C3）。
單位：L mm、Ec MPa、I mm⁴、w kN/m（= N/mm）、Pe N、e mm；撓度 mm（向下為正）。
關鍵：潛變同時作用於自重與預力（皆 ×(1+φ)）。預力反拱 w_eq = 8·Pe·a/L²。
"""
from dataclasses import dataclass
from .model import Section
from . import allowables


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
    d_LL_ok: bool     # δ_LL ≤ 該 case 之限值？
    d_LL_limit: float # 採用之撓度限值 mm
    defl_case: str    # 撓度限值情況（allowables.TW_DEFLECTION_DENOM 之鍵）


def deflection_analysis(L: float, Ec: float, section: Section,
                        w_DL: float, Pe: float, e: float, w_LL: float,
                        phi: float = 2.0, settlement: float = 5.0,
                        defl_case: str = "一般") -> DeflectionResult:
    """撓度與預拱。`defl_case` 依台灣 §7.3.12／§8.11.3／§9.1.7（2026-09-24 NLM 核）：
    '一般' L/800（預設，維持既有行為）／'市區人行' L/1000／'懸臂' L_c/300／'懸臂人行' L_c/375。
    ⚠ 懸臂情況請把 L 傳為**懸臂長度**。"""
    K = 5 * L**4 / (384 * Ec * section.I)         # mm/(N/mm)
    d_DL = K * w_DL
    w_eq = 8 * Pe * e / L**2                        # N/mm（= kN/m）
    d_PT = K * w_eq
    net_el = d_DL - d_PT
    LBR = d_PT / d_DL
    net_LT = net_el * (1 + phi)
    d_LL = K * w_LL
    camber = max(net_LT, 0.0) + settlement
    lim = allowables.deflection_limit_TW(L, defl_case)
    return DeflectionResult(K, d_DL, w_eq, d_PT, net_el, LBR, net_LT, d_LL,
                            camber, d_LL <= lim, lim, defl_case)
