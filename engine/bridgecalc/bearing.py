"""支承設計（E1）：疊層橡膠支承剪切應變、壓應力、上拔檢核。

對應公式卡_支承設計、算例_支承設計。位移量 Δ_S = 溫度+潛變+乾縮+預力縮短。
單位：力 kN、長度 mm、應力 MPa。
"""
from dataclasses import dataclass


@dataclass
class BearingResult:
    delta_s: float    # 設計位移量 mm
    gamma_s: float    # 剪切應變 Δ_S/h_rt
    sigma_TL: float   # 壓應力 MPa
    gamma_ok: bool    # γ_S ≤ 0.50
    no_uplift: bool   # R_min − R_LL > 0


def bearing_check(R_max_kN: float, R_min_kN: float, R_LL_kN: float,
                  delta_s: float, h_rt: float, plan_area_mm2: float,
                  gamma_limit: float = 0.50) -> BearingResult:
    """疊層橡膠支承關鍵檢核。delta_s 設計位移、h_rt 有效橡膠厚、plan_area 平面積。"""
    gamma_s = delta_s / h_rt
    sigma_TL = R_max_kN * 1e3 / plan_area_mm2          # MPa
    return BearingResult(delta_s, gamma_s, sigma_TL,
                         gamma_s <= gamma_limit, (R_min_kN - R_LL_kN) > 0)
