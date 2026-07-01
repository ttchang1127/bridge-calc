"""反解設計庫（階段 4）：把引擎從「驗算」升級為「設計」——閉式反解。

驗算 = 給定配置算應力/強度；設計 = 給定目標反解所需配置。
本模組收攏可符號導出的閉式反解（SymPy 為推導輔助，結果為純 Python，**零相依**）。
每式在 40m 參考橋上與正算閉環一致（見 tests / golden design_inverse）。
單位：力 N、長度 mm、應力 MPa（拉為正）、彎矩 kN·m；除非另註。
"""
import math

from .service import Pe_min_zero_tension   # 反解種子：底緣零拉所需最小 Pe（service.py）

__all__ = ["Pe_min_zero_tension", "min_tendon_groups",
           "required_drape", "min_section_modulus_Sb"]


def min_tendon_groups(Pe_req_N: float, strands_per: int, fpe_MPa: float,
                      Ap_strand: float = 140) -> int:
    """達所需有效預力 Pe 的最小鋼腱組數 = ⌈Pe / (股數·單股面積·f_pe)⌉。

    反解 Tendon 配置：Pe = n·股數·A_strand·f_pe → n = ⌈Pe/(股數·A_strand·f_pe)⌉。
    """
    per_group = strands_per * Ap_strand * fpe_MPa      # N/組
    return math.ceil(Pe_req_N / per_group)


def required_drape(LBR_target: float, w_DL: float, L: float, Pe_N: float) -> float:
    """達目標荷重平衡率 LBR 所需拋物線垂度 a（tendon_profile 反解）。

    由 w_eq = 8·Pe·a/L²、LBR = w_eq/w_DL 解得 **a = LBR·w_DL·L²/(8·Pe)**。
    w_DL[N/mm]、L[mm]、Pe[N] → a[mm]。
    """
    return LBR_target * w_DL * L ** 2 / (8 * Pe_N)


def min_section_modulus_Sb(Pe_N: float, A: float, e: float, M_kNm: float,
                           sigma_tension_limit: float = 0.0) -> float:
    """底緣不超拉所需最小底緣斷面模數 S_b（反解 service.stresses 的 σ_b）。

    σ_b = −Pe/A − Pe·e/S_b + M/S_b ≤ σ_lim（拉為正）
      → **S_b ≥ (M − Pe·e)/(σ_lim + Pe/A)**。M[kN·m]、σ_lim[MPa] → S_b[mm³]。
    與 Pe_min_zero_tension 為同一式 σ_b=0 的兩個反解（互為對偶）。
    """
    M = M_kNm * 1e6
    return (M - Pe_N * e) / (sigma_tension_limit + Pe_N / A)
