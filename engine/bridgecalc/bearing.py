"""支承設計（E1）：疊層橡膠支承完整檢核。

對應公式卡_支承設計、算例_支承設計。位移量 Δ_S = 溫度+潛變+乾縮+預力縮短。
檢核：剪切應變 γ_S、壓應力 σ_TL（≤ min(112 kgf/cm², 1.66GS)）、穩定性
（h_rt ≤ min(Lb,Wb)/3）、水平力 H_m（≤ R_min/5）、上拔。
單位：力 kN、長度 mm、應力 MPa；G 以 kgf/cm² 輸入（台灣 Method A）。
"""
from dataclasses import dataclass

KGF_CM2 = 0.0981   # kgf/cm² → MPa


@dataclass
class BearingResult:
    delta_s: float        # 設計位移量 mm
    gamma_s: float        # 剪切應變 Δ_S/h_rt
    sigma_TL: float       # 壓應力 MPa
    shape_S: float        # 形狀係數 S
    sigma_TL_limit: float # 壓應力上限 MPa = min(112 kgf, 1.66GS)
    H_m: float            # 水平力 kN
    gamma_ok: bool        # γ_S ≤ 0.50
    sigma_ok: bool        # σ_TL ≤ 上限
    stability_ok: bool    # h_rt ≤ min(Lb,Wb)/3
    H_ok: bool            # H_m ≤ R_min/5
    no_uplift: bool       # R_min − R_LL > 0


def bearing_check(R_max_kN: float, R_min_kN: float, R_LL_kN: float,
                  delta_s: float, h_rt: float, Lb: float, Wb: float,
                  te: float = 10.0, G_kgf: float = 8.0,
                  gamma_limit: float = 0.50) -> BearingResult:
    """疊層橡膠支承完整檢核。

    Lb×Wb 平面尺寸、te 單層橡膠厚、h_rt 有效橡膠總厚、G_kgf 剪力模數（kgf/cm²）。
    """
    A = Lb * Wb                                          # mm²
    gamma_s = delta_s / h_rt
    sigma_TL = R_max_kN * 1e3 / A                        # MPa
    S = (Lb * Wb) / (2 * te * (Lb + Wb))                 # 形狀係數
    limit_kgf = min(112.0, 1.66 * G_kgf * S)             # kgf/cm²
    sigma_TL_limit = limit_kgf * KGF_CM2                 # MPa
    G_MPa = G_kgf * KGF_CM2
    H_m = G_MPa * A * delta_s / h_rt / 1e3               # kN（= G·A·Δ/h_rt）
    stability_limit = min(Lb, Wb) / 3
    return BearingResult(
        delta_s, gamma_s, sigma_TL, S, sigma_TL_limit, H_m,
        gamma_s <= gamma_limit, sigma_TL <= sigma_TL_limit,
        h_rt <= stability_limit, H_m <= R_min_kN / 5,
        (R_min_kN - R_LL_kN) > 0)
