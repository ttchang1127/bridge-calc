"""推進工法施工設計（H4 ILM, Incremental Launching）：推進彎矩包絡、臨時置中預力、
滑動支承局部支壓、頂推力。

對應算例_推進工法施工設計、公式卡_推進工法施工設計。
ILM 特徵：梁段在後場澆置後逐段推進，每斷面在推進過程**交替承受正/負彎矩（包絡）**，
須以「臨時置中預力（e=0）」提供均勻壓力，使全斷面在最不利彎矩下仍保持受壓。
介面採工程單位：w[kN/m]、L[m]、M[kN·m]、A[mm²]、Z[mm³]、Pc/R[kN]、σ[MPa]。
"""
import math


def launching_cantilever_moment(w_kNpm: float, Lc_m: float) -> float:
    """推進中懸臂根部最大負彎矩 M⁻ = −w·Lc²/2（導梁未及次墩、最大懸臂時），kN·m。"""
    return -w_kNpm * Lc_m ** 2 / 2


def launching_span_moment(w_kNpm: float, L_m: float) -> float:
    """推進中跨中最大正彎矩 M⁺ ≈ w·L²/8，kN·m。"""
    return w_kNpm * L_m ** 2 / 8


def centric_prestress_required(M_pos_kNm: float, Zb_mm3: float, A_mm2: float,
                               sigma_res: float = 1.5) -> float:
    """臨時置中預力需求 Pc = (M⁺/Zb + σ_res)·A，kN。

    使底緣在最大正彎矩下仍保留 σ_res(MPa) 殘餘壓應力（不出現拉）。e=0 → 全斷面均勻壓。
    """
    f_demand = M_pos_kNm * 1e6 / Zb_mm3 + sigma_res    # MPa（底緣須抵銷的拉 + 殘餘壓）
    return f_demand * A_mm2 / 1e3                        # N → kN


def launching_bottom_stress(Pc_kN: float, A_mm2: float, M_pos_kNm: float,
                            Zb_mm3: float) -> float:
    """最大正彎矩工況底緣應力 σ_bot = −Pc/A + M⁺/Zb（MPa，壓為負）。"""
    return -Pc_kN * 1e3 / A_mm2 + M_pos_kNm * 1e6 / Zb_mm3


def n_tendons(Pc_kN: float, P_per_tendon_kN: float) -> int:
    """所需置中鋼腱束數 = ⌈Pc / 每束施拉力⌉。"""
    return math.ceil(Pc_kN / P_per_tendon_kN)


def jacking_force(mu_s: float, W_total_kN: float) -> float:
    """頂推力 F = μ_s·W_total（滑動摩擦；W_total = 梁體總重），kN。"""
    return mu_s * W_total_kN


def bearing_stress(R_kN: float, A_bearing_mm2: float) -> float:
    """滑動支承局部支壓應力 σ_ba = R / A_bearing（MPa）。"""
    return R_kN * 1e3 / A_bearing_mm2
