"""鋼腱線形設計（G1）：拋物線等效荷載法、荷重平衡率、摩擦損失、曲率半徑。

對應公式卡_鋼腱線形設計、算例_鋼腱線形設計。
平衡荷載法：拋物線鋼腱（垂度 a）對混凝土施加向上等效均佈荷重 w_eq = 8·P·a/L²，
用以抵銷靜載 w_DL；荷重平衡率 LBR = w_eq/w_DL 目標 0.75~0.90（施拉可暫時過平衡）。

單位：力 N、長度 mm（w_eq 回傳 N/mm，數值等於 kN/m）。
摩擦損失介面採規範慣用單位：k 以 /m、路徑長以 m。
"""
import math
from dataclasses import dataclass


def equivalent_load(P: float, a: float, L: float) -> float:
    """拋物線鋼腱等效均佈荷重 w_eq = 8·P·a/L²。P[N] a[mm] L[mm] → N/mm(=kN/m)。"""
    return 8 * P * a / L ** 2


def end_slope(a: float, L: float) -> float:
    """端部傾角 θ_end = 4a/L（對稱拋物線、兩端同高），rad。"""
    return 4 * a / L


def radius_of_curvature(a: float, L: float) -> float:
    """跨中最小曲率半徑 R = L²/(8a)，mm。"""
    return L ** 2 / (8 * a)


def balance_ratio(w_eq: float, w_DL: float) -> float:
    """荷重平衡率 LBR = w_eq / w_DL。"""
    return w_eq / w_DL


def friction_loss(alpha: float, mu: float, k_per_m: float, L_path_m: float) -> float:
    """後張摩擦損失率 ΔP/P = 1 − e^−(μ·α + k·L)。α 累積角變化[rad]、k[/m]、L[m]。"""
    return 1 - math.exp(-(mu * alpha + k_per_m * L_path_m))


@dataclass
class TendonProfileResult:
    a: float               # 垂度 mm
    theta_end: float       # 端部傾角 rad
    R: float               # 曲率半徑 mm
    w_eq_transfer: float   # 施拉時等效荷重 N/mm(=kN/m)
    w_eq_service: float    # 服務時等效荷重 N/mm(=kN/m)
    LBR_transfer: float    # 施拉荷重平衡率
    LBR_service: float     # 服務荷重平衡率
    fric_single_end: float # 單端張拉、遠端（全長）損失率
    fric_dual_mid: float   # 雙端張拉、跨中損失率
    R_ok: bool             # R ≥ R_min


def tendon_profile(Pi: float, Pe: float, a: float, L: float, w_DL: float,
                   mu: float = 0.25, k_per_m: float = 0.003,
                   R_min: float = 3000) -> TendonProfileResult:
    """完整 G1 線形驗算。

    Pi/Pe[N]：施拉/有效預力；a[mm]：跨中垂度（= e_m − e_support，端部偏心 0 時即 e_m）；
    L[mm]：跨徑；w_DL[N/mm]：靜載（簡支可由 8(M_DC+M_DW)/L² 反推）。
    摩擦：單端張拉取全長 α=2θ_end、路徑 L；雙端張拉取跨中 α=θ_end、路徑 L/2。
    """
    th = end_slope(a, L)
    R = radius_of_curvature(a, L)
    w_t = equivalent_load(Pi, a, L)
    w_s = equivalent_load(Pe, a, L)
    L_m = L / 1000.0
    return TendonProfileResult(
        a=a, theta_end=th, R=R,
        w_eq_transfer=w_t, w_eq_service=w_s,
        LBR_transfer=balance_ratio(w_t, w_DL),
        LBR_service=balance_ratio(w_s, w_DL),
        fric_single_end=friction_loss(2 * th, mu, k_per_m, L_m),
        fric_dual_mid=friction_loss(th, mu, k_per_m, L_m / 2.0),
        R_ok=(R >= R_min),
    )
