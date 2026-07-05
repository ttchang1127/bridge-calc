"""耐震設計閉式（台灣公路橋梁耐震設計規範，單軌）。

把耐震線四張公式卡的**閉式**收攏為可回歸的純 Python（零相依）：
  S1 落橋防止：min N_L（式8-10）、地盤相對變位 u_G、防落長度需求 L_N
  S2 隔震消能：LRB 雙線性 + 等效線性化**迭代**（式7-1~7-5、表3-1 阻尼修正、式C7）
  S3 橋墩韌性：超強彎矩 M_p、容量剪力 V_u、圍束 ρ_s（式5-5/5-6）/A_sh（式5-7/5-8）、間距/塑鉸長度
  S5 液狀化：土壤參數折減係數 D_E 查表（表8-1）

> 與箱梁主引擎（雙軌 HS20/HL-93）不同，耐震規範屬**台灣單軌**；本模組
  沿用各條文原生單位（詳見各函式 docstring），不做跨軌換算。
> 每式在對應算例上閉環一致（見 golden seismic_S1/S2/S3/S5）。
"""
import math
from dataclasses import dataclass

__all__ = [
    # S1 落橋防止
    "min_falloff_length", "ground_relative_displacement", "required_falloff_length",
    "restrainer_yield_strength", "EPSILON_G",
    # S2 隔震
    "isolation_bilinear", "damping_correction_B1", "damping_correction_BS",
    "isolation_design", "IsolationResult",
    # S3 橋墩韌性
    "overstrength_moment", "capacity_shear", "rho_s_circular", "Ash_rectangular",
    "confinement_spacing_limit", "plastic_hinge_length",
    # S5 液狀化
    "liquefaction_reduction_DE",
]

# ────────────────────────────────────────────────────────────────
# S1 落橋防止系統（§8.5）  單位：長度 cm、L/H m、S 度、力 kN
# ────────────────────────────────────────────────────────────────

# 地盤變位係數 ε_G（表 §8.5，依地盤類別）
EPSILON_G = {"第一類": 0.0025, "第二類": 0.00375, "第三類": 0.005, "臺北盆地": 0.00625}


def min_falloff_length(L: float, H: float, S: float = 0.0) -> float:
    """最小梁端防落長度 min N_L（式8-10），單位 cm。

    min N_L = (50 + 0.25·L + 1.0·H)·(1 + S²/8000)
    L=跨徑[m]、H=基面起算下部結構高[m]、S=橋墩斜角[度]。
    """
    return (50 + 0.25 * L + 1.0 * H) * (1 + S ** 2 / 8000)


def ground_relative_displacement(soil_class: str, L_e: float,
                                 ratio_SIII_SII: float) -> float:
    """地震引致相鄰墩間地盤水平相對變位 u_G，單位 cm。

    u_G = ε_G · L_e · (S_III,S / S_II,S)
    soil_class ∈ EPSILON_G；L_e=影響防落長度之下部結構間距[cm]；
    ratio=等級Ⅲ/等級Ⅱ短週期譜加速度係數比。
    """
    return EPSILON_G[soil_class] * L_e * ratio_SIII_SII


def required_falloff_length(min_NL: float, u_R: float = 0.0, u_G: float = 0.0,
                            movable: bool = True) -> float:
    """防落長度需求 L_N，單位 cm。

    活動支承：L_N ≥ max(min N_L, u_R + u_G)（式8-10 兩條件取大）
    固定支承（movable=False）：L_N ≥ min N_L。
    u_R=等級Ⅲ梁墩相對變位[cm]、u_G=地盤相對變位[cm]。
    """
    if not movable:
        return min_NL
    return max(min_NL, u_R + u_G)


def restrainer_yield_strength(R_d: float) -> float:
    """落橋防止裝置設計降伏強度下限 F_y = 1.5·R_d（§8.5）。R_d=支承靜載反力[kN]→F_y[kN]。"""
    return 1.5 * R_d


# ────────────────────────────────────────────────────────────────
# S2 隔震與消能設計（第7章）  單位：力 kN、長度 m、勁度 kN/m、週期 s
# ────────────────────────────────────────────────────────────────

def isolation_bilinear(Q_d: float, K_d: float, D_d: float):
    """LRB 雙線性：有效（割線）勁度與等效阻尼比。

    K_eff = Q_d/D_d + K_d（式）；E_D ≈ 4·Q_d·D_d（遲滯迴圈，D_y 略）；
    ξ_eq = E_D/(2π·K_eff·D_d²) = (2/π)·Q_d/(K_eff·D_d)（式7-4）。
    Q_d[kN]、K_d[kN/m]、D_d[m] → (K_eff[kN/m], ξ_eq[-])。
    """
    K_eff = Q_d / D_d + K_d
    E_D = 4 * Q_d * D_d
    xi_eq = E_D / (2 * math.pi * K_eff * D_d ** 2)
    return K_eff, xi_eq


def damping_correction_B1(xi: float) -> float:
    """一秒週期阻尼修正係數 B_1（表3-1，線性內插）。xi 為分數（0.15=15%）。

    控制點(ξ%,B1)：(<2,0.80)(5,1.00)(10,1.25)(≥20,1.50)。
    """
    p = xi * 100.0
    if p <= 2:
        return 0.80
    if p <= 5:
        return 0.80 + (1.00 - 0.80) * (p - 2) / (5 - 2)
    if p <= 10:
        return 1.00 + (1.25 - 1.00) * (p - 5) / (10 - 5)
    if p <= 20:
        return 1.25 + (1.50 - 1.25) * (p - 10) / (20 - 10)
    return 1.50


def damping_correction_BS(xi: float) -> float:
    """短週期阻尼修正係數 B_S（表3-1，線性內插）。控制點(ξ%,BS)：(<2,0.80)(5,1.00)(10,1.33)(≥20,1.60)。"""
    p = xi * 100.0
    if p <= 2:
        return 0.80
    if p <= 5:
        return 0.80 + (1.00 - 0.80) * (p - 2) / (5 - 2)
    if p <= 10:
        return 1.00 + (1.33 - 1.00) * (p - 5) / (10 - 5)
    if p <= 20:
        return 1.33 + (1.60 - 1.33) * (p - 10) / (20 - 10)
    return 1.60


@dataclass
class IsolationResult:
    """隔震等效線性化迭代結果。"""
    D_d: float          # 隔震元件設計位移 [m]
    T_e: float          # 系統有效週期 [s]
    K_eff: float        # 隔震元件有效勁度 [kN/m]
    xi_e: float         # 系統等效阻尼比 [-]
    B1: float           # 阻尼修正係數（表3-1）
    S_a: float          # 譜加速度係數 S_a,II [-]
    V_b_secant: float   # 設計剪力 K_eff·D_d（式C7-3）[kN]
    V_b_bilinear: float # 設計剪力 Q_d+K_d·D_d（式C7-4）[kN]
    iterations: int     # 收斂迭代次數


def isolation_design(W: float, Q_d: float, K_d: float, S_II_1: float,
                     K_p: float = None, T0_II: float = None, S_II_S: float = None,
                     taipei_basin: bool = False, g: float = 9.81,
                     tol: float = 1e-4, max_iter: int = 100,
                     s_D0: float = 0.20) -> IsolationResult:
    """隔震橋等效線性化**迭代**（解說 C7.3）。

    因 K_eff、ξ_eq 皆為設計位移 D_d 的函數，須迭代收斂：
      s_D → D_d → K_eff/ξ_eq → (串聯 K_e) → T_e → B_1(ξ_e) → S_a,II → 新 s_D
    收斂後 D_d、T_e、V_b。W[kN]、Q_d[kN]、K_d[kN/m]、K_p=下部結構勁度[kN/m]（None=剛性墩）。

    S_a,II（T_e>T0_II，1秒段）：一般 S_II,1/(B_1·T_e)（式7-1b）；
                               臺北盆地 T0_II·S_II,S/(B_1·T_e)（式7-1c）。
    剛性墩：D_d=s_D、K_e=K_eff、ξ_e=ξ_eq。柔性墩：D_P 由式C7-2、K_e 串聯式7-3、
    ξ_e 保守取 ξ_eq（忽略下部/基礎耗能，見公式卡 S2 §五註）。
    """
    s_D = s_D0
    it = 0
    K_eff = xi_e = T_e = B1 = S_a = D_d = 0.0
    for it in range(1, max_iter + 1):
        if K_p is None:
            D_d = s_D
            K_eff, xi_eq = isolation_bilinear(Q_d, K_d, D_d)
            K_e = K_eff
            xi_e = xi_eq
        else:
            D_P = (Q_d + K_d * s_D) / (K_p + K_d)   # 墩頂位移（式C7-2）
            D_d = s_D - D_P                          # 隔震器位移（式C7-1）
            K_eff, xi_eq = isolation_bilinear(Q_d, K_d, D_d)
            K_e = K_eff * K_p / (K_eff + K_p)        # 串聯有效勁度（式7-3）
            xi_e = xi_eq                             # 保守（忽略下部/基礎耗能）
        T_e = 2 * math.pi * math.sqrt(W / (g * K_e))
        B1 = damping_correction_B1(xi_e)
        if taipei_basin:
            S_a = T0_II * S_II_S / (B1 * T_e)        # 式7-1c
        else:
            S_a = S_II_1 / (B1 * T_e)                # 式7-1b
        new_sD = S_a * T_e ** 2 / (4 * math.pi ** 2) * g   # 式7-1a
        if abs(new_sD - s_D) < tol:
            s_D = new_sD
            break
        s_D = new_sD
    if K_p is None:
        D_d = s_D
    else:
        D_d = s_D - (Q_d + K_d * s_D) / (K_p + K_d)
    return IsolationResult(
        D_d=D_d, T_e=T_e, K_eff=K_eff, xi_e=xi_e, B1=B1, S_a=S_a,
        V_b_secant=K_eff * D_d, V_b_bilinear=Q_d + K_d * D_d, iterations=it)


# ────────────────────────────────────────────────────────────────
# S3 橋墩韌性耐震設計（§4.2 / 5.3）  單位：規範原生 kgf/cm²、cm、tf·m
# ────────────────────────────────────────────────────────────────

def overstrength_moment(M_n: float, overstrength: float = 1.3) -> float:
    """最大可能（超強）彎矩 M_p = 1.3·M_n（§4.2.1）。M_n 任意彎矩單位，回同單位。"""
    return overstrength * M_n


def capacity_shear(sum_Mp: float, L_c: float) -> float:
    """容量設計剪力 V_u = ΣM_p / L_c。單柱懸臂 sum_Mp=M_p、L_c=H；雙鉸取兩端 M_p 和。

    M_p[tf·m]、L_c[m] → V_u[tf]（單位需自洽）。
    """
    return sum_Mp / L_c


def rho_s_circular(fc: float, fyh: float, Ag: float, Ac: float, Pe: float) -> float:
    """圓柱塑鉸區螺箍體積比 ρ_s = max(式5-5幾何式, 式5-6軸力式)。

    式5-5：0.45·(f'c/fyh)·(Ag/Ac − 1)
    式5-6：0.12·(f'c/fyh)·(0.5 + 1.25·Pe/(f'c·Ag))
    f'c,fyh[kgf/cm²]、Ag,Ac[cm²]、Pe[kgf] → ρ_s[-]（大軸力常由式5-6控制）。
    """
    r55 = 0.45 * (fc / fyh) * (Ag / Ac - 1)
    r56 = 0.12 * (fc / fyh) * (0.5 + 1.25 * Pe / (fc * Ag))
    return max(r55, r56)


def Ash_rectangular(a: float, hc: float, fc: float, fyh: float,
                    Ag: float, Ac: float, Pe: float) -> float:
    """矩柱塑鉸區單方向橫箍總斷面積 A_sh = max(式5-7, 式5-8)，單位 cm²（兩主軸各算）。

    式5-7：0.30·a·hc·(f'c/fyh)·(Ag/Ac − 1)
    式5-8：0.12·a·hc·(f'c/fyh)·(0.5 + 1.25·Pe/(f'c·Ag))
    a=箍筋垂直間距[cm]、hc=柱心所考慮方向尺寸[cm]。
    """
    a57 = 0.30 * a * hc * (fc / fyh) * (Ag / Ac - 1)
    a58 = 0.12 * a * hc * (fc / fyh) * (0.5 + 1.25 * Pe / (fc * Ag))
    return max(a57, a58)


def confinement_spacing_limit(short_side: float, d_b: float) -> float:
    """圍束筋垂直間距上限 a ≤ min(15cm, 柱短邊/4, 6·d_b)（§5.3.5第3款）。單位 cm。"""
    return min(15.0, short_side / 4.0, 6.0 * d_b)


def plastic_hinge_length(col_depth: float, L_c: float, floor: float = 45.0) -> float:
    """塑鉸區配置長度 ℓ₀ = max(沿剪力方向柱深, ℓ_c/6, 45cm)（§5.3.5第1款）。單位 cm。"""
    return max(col_depth, L_c / 6.0, floor)


# ────────────────────────────────────────────────────────────────
# S5 液狀化（§8.1）  無因次折減係數查表
# ────────────────────────────────────────────────────────────────

def liquefaction_reduction_DE(F_L: float, x: float, R_s: float) -> float:
    """液化土層土壤參數折減係數 D_E（表8-1）。

    依 液化抵抗率 F_L、距地表深度 x[m]、抵抗液化剪應力強度比 R_s 查表。
    F_L≥1.0 → 不折減(1.0)；x>20m 超出表範圍亦回 1.0（不折減）。
    回傳側向地盤反力係數與極限反力之折減乘數（0=參數設零）。
    """
    if F_L >= 1.0 or x > 20:
        return 1.0
    shallow = x <= 10                      # 0≤x≤10 為淺層
    dense = R_s > 0.3                      # 密實砂折減較輕
    if F_L <= 1.0 / 3.0:                   # 第一級
        if shallow:
            return 1.0 / 6.0 if dense else 0.0
        return 1.0 / 3.0
    if F_L <= 2.0 / 3.0:                   # 第二級
        if shallow:
            return 2.0 / 3.0 if dense else 1.0 / 3.0
        return 2.0 / 3.0
    # 第三級 2/3 < F_L ≤ 1.0
    if shallow:
        return 1.0 if dense else 2.0 / 3.0
    return 1.0
