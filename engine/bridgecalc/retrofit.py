"""補強設計閉式（中國公路橋梁加固設計規範 JTG/T J22，單軌）。

把廣義補強線的正斷面抗彎**承載力式**收攏為純 Python（零相依）：
  R1 碳纖維CFRP：允許拉應變 [ε_f]（剝離控制 κ_m）、案② 承載力（式6-42）
  R2 外貼鋼板：鋼板應變 ε_sp（二次受力）、承載力（式6-26）、粘貼延伸長度 l_p（式6-37）
  R4 增大截面：新增筋應變 ε_s2、承載力（式6-2/6-3）
  共用：開裂換算斷面 x₁/I_cr、一期初始應變 ε_c1（二次受力核心）

> 與箱梁主引擎（雙軌 HS20/HL-93）不同，加固規範屬**中國 JTG 單軌**；沿用原文
  符號與設計值。每式在對應算例上閉環一致（見 golden retrofit_R1/R2/R4）。
> 單位：應力 MPa、長度 mm、力 N；彎矩回傳 kN·m。壓為正（承載力式取正號體系）。
"""
import math
from dataclasses import dataclass

__all__ = [
    "cracked_na_depth", "cracked_inertia", "initial_concrete_strain",
    "cfrp_km1", "cfrp_allowable_strain", "xi_fb",
    "cfrp_moment_capacity", "CFRPResult",
    "plate_moment_capacity", "plate_dev_length", "PlateResult",
    "enlargement_moment_capacity", "EnlargeResult",
]


# ──────────────────── 共用：開裂換算斷面（二次受力）────────────────────

def cracked_na_depth(b: float, h0: float, As: float, alpha_Es: float,
                     As_comp: float = 0.0, as_comp: float = 0.0) -> float:
    """加固前開裂換算斷面受壓區高 x₁ = √(A₁²+B₁) − A₁（彈性，矩形）。

    A₁=α_Es(A_s+A′_s)/b、B₁=2α_Es(A_s·h₀+A′_s·a′_s)/b。α_Es=E_s/E_c。
    """
    A1 = alpha_Es * (As + As_comp) / b
    B1 = 2 * alpha_Es * (As * h0 + As_comp * as_comp) / b
    return math.sqrt(A1 ** 2 + B1) - A1


def cracked_inertia(b: float, x1: float, As: float, h0: float, alpha_Es: float,
                    As_comp: float = 0.0, as_comp: float = 0.0) -> float:
    """加固前開裂換算斷面慣性矩 I_cr = b·x₁³/3 + α_Es·A_s(h₀−x₁)² + α_Es·A′_s(x₁−a′_s)²。"""
    return (b * x1 ** 3 / 3 + alpha_Es * As * (h0 - x1) ** 2
            + alpha_Es * As_comp * (x1 - as_comp) ** 2)


def initial_concrete_strain(Md1_kNm: float, x1: float, Ec: float, Icr: float) -> float:
    """一期荷載下混凝土初始應變 ε_c1 = M_d1·x₁/(E_c·I_cr)。M_d1[kN·m]、x₁/I_cr[mm]。"""
    return Md1_kNm * 1e6 * x1 / (Ec * Icr)


# ──────────────────── R1 碳纖維 CFRP（式6-42/6-44）────────────────────

def cfrp_km1(n_f: int, E_f: float, t_f: float) -> float:
    """CFRP 層數/厚度折減 κ_m1（式6-44）。n_f 層、E_f[MPa]、t_f[mm]。"""
    v = n_f * E_f * t_f
    return 1 - v / 428000 if v <= 214000 else 1070000 / v


def cfrp_allowable_strain(n_f: int, E_f: float, t_f: float, eps_fu: float,
                          km2: float = 0.85) -> float:
    """CFRP 允許拉應變 [ε_f] = min(κ_m·ε_fu, ⅔ε_fu, 0.007)。κ_m=min(κ_m1,κ_m2)≤0.9。

    km2 環境折減（碳纖維 0.85）。剝離控制，防高估承載。
    """
    km = min(cfrp_km1(n_f, E_f, t_f), km2, 0.9)
    return min(km * eps_fu, 2 / 3 * eps_fu, 0.007)


def xi_fb(eps_cu: float, eps_f: float, eps_1: float) -> float:
    """界限相對受壓區高 ξ_fb = 0.8·ε_cu/(ε_cu+[ε_f]+ε₁)（CFRP 達限與混凝土壓壞同時）。"""
    return 0.8 * eps_cu / (eps_cu + eps_f + eps_1)


@dataclass
class CFRPResult:
    eps_f_allow: float   # 允許拉應變 [ε_f]
    xi_fb: float         # 界限相對受壓區高
    x: float             # 受壓區高（案②平衡）mm
    Mu_kNm: float        # 補強後抗彎承載力 kN·m
    case2: bool          # x ≤ ξ_fb·h → 案②（CFRP 達 [ε_f]）成立


def cfrp_moment_capacity(b: float, h: float, h0: float, fcd: float, eps_cu: float,
                         As: float, fsd: float, Af: float, Ef: float,
                         eps_f_allow: float, eps_1: float) -> CFRPResult:
    """CFRP 補強後正斷面抗彎承載力（案②，式6-42；單筋略受壓筋）。

    案②（x ≤ ξ_fb·h，CFRP 達 [ε_f]）：
      M_u = f_sd·A_s(h₀−0.5ξ_fb·h) + E_f·[ε_f]·A_f·h(1−0.5ξ_fb)
    x 由軸力平衡 f_cd·b·x = f_sd·A_s + E_f·[ε_f]·A_f 判別案別。h0=鋼筋深、h=斷面高(CFRP在底)。
    """
    xfb = xi_fb(eps_cu, eps_f_allow, eps_1)
    x = (fsd * As + Ef * eps_f_allow * Af) / (fcd * b)
    Mu = fsd * As * (h0 - 0.5 * xfb * h) + Ef * eps_f_allow * Af * h * (1 - 0.5 * xfb)
    return CFRPResult(eps_f_allow, xfb, x, Mu / 1e6, x <= xfb * h)


# ──────────────────── R2 外貼鋼板（式6-26/6-28/6-35/6-37）────────────────────

@dataclass
class PlateResult:
    x: float             # 受壓區高 mm
    Mu_kNm: float        # 補強後抗彎承載力 kN·m
    plate_yields: bool   # 鋼板是否降伏（σ_sp 需求 ≥ f_sp）
    sigma_sp: float      # 鋼板工作應力 MPa（≤ f_sp）


def plate_moment_capacity(b: float, h: float, h0: float, fcd1: float, eps_cu: float,
                          As: float, fsd: float, Asp: float, fsp: float, Esp: float,
                          x1: float, eps_c1: float, beta: float = 0.8,
                          a_s: float = None) -> PlateResult:
    """外貼鋼板補強後正斷面抗彎承載力（式6-26；單筋略受壓筋）。

    鋼板應變 ε_sp = ε_cu(βh−x)/x − ε_c1(h−x₁)/x₁（式6-35，二次受力），σ_sp=E_sp·ε_sp≤f_sp；
    軸力平衡 f_cd1·b·x = f_sd·A_s + σ_sp·A_sp（式6-28，σ_sp/x 耦合，迭代）；
    M_u = f_cd1·b·x(h₀−x/2) + σ_sp·A_sp·a_s（式6-26，a_s=鋼筋至鋼板距=h−h₀）。
    """
    if a_s is None:
        a_s = h - h0
    sigma_sp = fsp
    x = eps_sp = 0.0
    for _ in range(50):
        x = (fsd * As + sigma_sp * Asp) / (fcd1 * b)
        eps_sp = eps_cu * (beta * h - x) / x - eps_c1 * (h - x1) / x1
        new = min(Esp * eps_sp, fsp)
        if abs(new - sigma_sp) < 0.01:
            sigma_sp = new
            break
        sigma_sp = new
    Mu = fcd1 * b * x * (h0 - x / 2) + sigma_sp * Asp * a_s
    return PlateResult(x, Mu / 1e6, Esp * eps_sp >= fsp, sigma_sp)


def plate_dev_length(fsp: float, Asp: float, tau_p: float, b_p: float) -> float:
    """粘貼延伸長度 l_p = f_sp·A_sp/(τ_p·b_p) + 300（式6-37），mm。"""
    return fsp * Asp / (tau_p * b_p) + 300


# ──────────────────── R4 增大截面（式6-2/6-3/6-10）────────────────────

@dataclass
class EnlargeResult:
    x: float               # 受壓區高 mm
    Mu_kNm: float          # 補強後抗彎承載力 kN·m
    added_bar_yields: bool # 新增鋼筋是否降伏
    sigma_s2: float        # 新增鋼筋工作應力 MPa（≤ f_sd2）


def enlargement_moment_capacity(b2: float, h02: float, h0: float, fcd1: float,
                                eps_cu: float, As1: float, fsd1: float,
                                As2: float, fsd2: float, Es2: float,
                                x1: float, eps_c1: float, beta: float = 0.8,
                                As_comp: float = 0.0, fsd_comp: float = 0.0,
                                as_comp: float = 0.0) -> EnlargeResult:
    """增大截面補強後正斷面抗彎承載力（式6-2）。

    新增筋應變 ε_s2 = ε_cu(β·h₀₂−x)/x − ε_c1(h₀₂−x₁)/x₁（式6-10，應變滯後），σ_s2=E_s2·ε_s2≤f_sd2；
    軸力平衡 f_cd1·b₂·x = f_sd1·A_s1 − f′_sd·A′_s + σ_s2·A_s2（式6-3，迭代）；
    M_u = f_cd1·b₂·x(h₀−x/2) + f′_sd·A′_s(h₀−a′_s)（式6-2，h₀=A_s1+A_s2 合力點）。
    """
    sigma_s2 = fsd2
    x = eps_s2 = 0.0
    for _ in range(50):
        x = (fsd1 * As1 - fsd_comp * As_comp + sigma_s2 * As2) / (fcd1 * b2)
        eps_s2 = eps_cu * (beta * h02 - x) / x - eps_c1 * (h02 - x1) / x1
        new = min(Es2 * eps_s2, fsd2)
        if abs(new - sigma_s2) < 0.01:
            sigma_s2 = new
            break
        sigma_s2 = new
    Mu = fcd1 * b2 * x * (h0 - x / 2) + fsd_comp * As_comp * (h0 - as_comp)
    return EnlargeResult(x, Mu / 1e6, Es2 * eps_s2 >= fsd2, sigma_s2)
