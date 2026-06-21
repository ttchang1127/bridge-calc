"""後張預力損失（台灣／AASHTO 集總法）。對應公式卡 B1（短期）、B2（長期）。

★ 本模組捕捉「增配的非線性耦合」——手動傳播時最易算錯的一環：
   f_cgp = Pi(1/A + e²/I) − M_D·e/I
   其中自重項 M_D·e/I 為**常數、不隨股數放大**，故增配（股數↑→Pi↑）會使
   f_cgp 升幅 > 股數升幅 → ES 與潛變損失增加 → fpe 下降，Pe 非線性。
   （見 算例_後張箱梁服務性應力驗算 §十觀察 3）
"""
from dataclasses import dataclass
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
                   friction: float = 117.0, shrink: float = 32.0,
                   relax: float = 10.0, Ep_Eci: float = 7.33) -> LossResult:
    """以集總法計算總損失與有效預力 Pe。

    friction/shrink/relax 為幾何/環境相關（與股數無關），預設取 40m 參考橋值。
    ES 與 creep 由 f_cgp 推算 → 隨配置自動變動（非線性耦合的來源）。
    """
    Pi = tendon.Pi
    e = tendon.e
    M_D = M_D_kNm * 1e6
    M_SDL = M_SDL_kNm * 1e6

    fcgp = Pi / section.A + Pi * e**2 / section.I - M_D * e / section.I  # MPa
    fcir = fcgp
    fcds = M_SDL * e / section.I

    N = tendon.n_tendons
    ES = (N - 1) / (2 * N) * Ep_Eci * fcgp          # 後張順序張拉
    short = friction + ES                            # 錨具滑移跨中=0（影響長度未達）

    creep = 12.0 * fcir - 7.0 * fcds                 # 台灣法（MPa）
    long = creep + shrink + relax

    total = short + long
    fpe = tendon.fpj - total
    Pe = fpe * tendon.Aps                            # N
    return LossResult(fcgp, ES, creep, friction, shrink, relax,
                      short, long, total, total / tendon.fpj, fpe, Pe)
