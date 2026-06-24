"""錨碇區爆裂力（F1）：整體區爆裂力 Tburst、爆裂筋、剝落力。

對應公式卡_錨碇區爆裂力設計、算例_錨碇區爆裂力設計。逐束法（等間距佈置）。
單位：力 kN、長度 mm、應力 MPa。
"""
from dataclasses import dataclass


@dataclass
class AnchorageResult:
    Pu: float          # 每束設計力 kN
    Tburst: float      # 每束爆裂力 kN
    sum_Tburst: float  # 單腹版群爆裂力 kN
    As_burst: float    # 爆裂筋 mm²（每束）
    Fspall: float      # 剝落力 kN（腹版群）
    As_spall: float    # 剝落筋 mm²


def anchorage_check(Pi_total_kN: float, n_tendons: int, a_plate: float,
                    h_diaph: float, n_per_web: int, fy: float = 400.0,
                    load_factor: float = 1.2) -> AnchorageResult:
    """端橫隔版錨碇區（逐束法）。a_plate 錨板高、h_diaph 端橫隔版全高。"""
    Pu = load_factor * Pi_total_kN / n_tendons
    Tburst = 0.25 * Pu * (1 - a_plate / h_diaph)
    sum_Tb = n_per_web * Tburst
    As_burst = 1.2 * Tburst * 1e3 / (0.85 * fy)
    sum_Pu = n_per_web * Pu
    Fspall = 0.02 * sum_Pu
    As_spall = Fspall * 1e3 / (0.6 * fy)
    return AnchorageResult(Pu, Tburst, sum_Tb, As_burst, Fspall, As_spall)


def spiral_local_bearing(Pu_kN: float, Pult0_kN: float, flat: float,
                         A_core: float, s: float, Dsp: float):
    """局部錨碇區螺旋圍束承壓強度（STEP 1，AASHTO 局部區）。

    P_ult = P_ult0 + 4.1·f_lat·A_core·(1−s/Dsp)²；f_lat 圍束壓力 MPa、A_core 核心面積 mm²、
    s 螺距、Dsp 螺圈外徑。回傳 (P_ult kN, 承壓餘裕, ok)。
    """
    Pult = Pult0_kN + 4.1 * flat * A_core * (1 - s / Dsp) ** 2 / 1e3
    return Pult, Pult / Pu_kN, Pult >= Pu_kN
