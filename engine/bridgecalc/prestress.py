"""後張預力損失（台灣／AASHTO 集總法）。對應公式卡 B1（短期）、B2（長期）。

★ 本模組捕捉「增配的非線性耦合」——手動傳播時最易算錯的一環：
   f_cgp = Pi(1/A + e²/I) − M_D·e/I
   其中自重項 M_D·e/I 為**常數、不隨股數放大**，故增配（股數↑→Pi↑）會使
   f_cgp 升幅 > 股數升幅 → ES 與潛變損失增加 → fpe 下降，Pe 非線性。
   （見 算例_後張箱梁服務性應力驗算 §十觀察 3）
"""
from dataclasses import dataclass
from math import exp
from .model import Section, Tendon


@dataclass
class LossResult:
    fcgp: float       # 鋼腱形心處混凝土應力 MPa
    ES: float         # 彈性縮短損失 MPa
    creep: float      # 潛變損失 MPa
    friction: float   # 摩擦損失 MPa
    shrink: float     # 乾縮損失 MPa
    relax: float      # 鬆弛損失 MPa
    short: float      # 短期合計 MPa
    long: float       # 長期合計 MPa
    total: float      # 總損失 MPa
    loss_pct: float   # 總損失率（/fpj）
    fpe: float        # 有效預力應力 MPa
    Pe: float         # 有效預力 N


def compute_losses(tendon: Tendon, section: Section,
                   M_D_kNm: float, M_SDL_kNm: float,
                   mu: float = 0.25, K: float = 0.003, alpha: float = 0.111,
                   x_ctrl: float = 20.0, RH: float = 75.0,
                   relax: float = 10.0, Ep_Eci: float = 7.33) -> LossResult:
    """計算總損失與有效預力 Pe（各損失由第一原理參數推導）。

    摩擦  ΔfpF = fpj·(1−e^(−(K·x+μ·α)))           （公式卡_後張預力短期損失 §一）
    乾縮  SH  = 0.8·(1195−10.55·RH)（kgf/cm²）×0.0981 （長期損失 §二，台灣式）
    ES、creep 由 f_cgp 推算 → 隨配置自動變動（非線性耦合的來源）。
    鬆弛  relax：低鬆弛型 ~8–10 MPa（AASHTO 5.9.5.4.2c 簡化 8；本庫取 10），
                量小且公式分歧，保留為參數。錨具滑移跨中=0（影響長度未達 L/2）。
    參數 μ/K/α/x_ctrl/RH 預設為 40m 參考橋值（雙端張拉、跨中控制 x=L/2=20m）。
    """
    Pi = tendon.Pi
    e = tendon.e
    M_D = M_D_kNm * 1e6
    M_SDL = M_SDL_kNm * 1e6

    fcgp = Pi / section.A + Pi * e**2 / section.I - M_D * e / section.I  # MPa
    fcir = fcgp
    fcds = M_SDL * e / section.I

    friction = tendon.fpj * (1 - exp(-(K * x_ctrl + mu * alpha)))   # ΔfpF（摩擦+偏折）
    shrink = 0.8 * (1195 - 10.55 * RH) * 0.0981                     # SH（kgf/cm²→MPa）

    N = tendon.n_tendons
    ES = (N - 1) / (2 * N) * Ep_Eci * fcgp          # 後張順序張拉
    short = friction + ES                            # 錨具滑移跨中=0（影響長度未達）

    creep = 12.0 * fcir - 7.0 * fcds                 # 台灣法（MPa，12/7 為無因次係數）
    long = creep + shrink + relax

    total = short + long
    fpe = tendon.fpj - total
    Pe = fpe * tendon.Aps                            # N
    return LossResult(fcgp, ES, creep, friction, shrink, relax,
                      short, long, total, total / tendon.fpj, fpe, Pe)
