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
                   relax: float = 10.0, Ep_Eci: float = 7.33,
                   fric_ratio: float = None) -> LossResult:
    """計算總損失與有效預力 Pe（各損失由第一原理參數推導）。

    摩擦  ΔfpF = fpj·(1−e^(−(K·x+μ·α)))           （公式卡_後張預力短期損失 §一）
    乾縮  SH  = 0.8·(1195−10.55·RH)（kgf/cm²）×0.0981 （長期損失 §二，台灣式）
    ES、creep 由 f_cgp 推算 → 隨配置自動變動（非線性耦合的來源）。
    鬆弛  relax：低鬆弛型 ~8–10 MPa（AASHTO 5.9.5.4.2c 簡化 8；本庫取 10），
                量小且公式分歧，保留為參數。錨具滑移跨中=0（影響長度未達 L/2）。
    參數 μ/K/α/x_ctrl/RH 預設為 40m 參考橋值（雙端張拉、跨中控制 x=L/2=20m）。
    fric_ratio：直接指定摩擦損失率（見 tendon_profile.friction_at，依張拉端配置與
    線形算出該控制斷面的值）；不給則沿用 α/x_ctrl 的單點式。
    """
    Pi = tendon.Pi
    e = tendon.e
    M_D = M_D_kNm * 1e6
    M_SDL = M_SDL_kNm * 1e6

    fcgp = Pi / section.A + Pi * e**2 / section.I - M_D * e / section.I  # MPa
    fcir = fcgp
    fcds = M_SDL * e / section.I

    # fric_ratio 給定時直接採用（由 tendon_profile.friction_at 依張拉端配置與線形算得），
    # 否則走原本的單點 α/x_ctrl 式——不傳即與既有結果完全相同。
    fr = (1 - exp(-(K * x_ctrl + mu * alpha))) if fric_ratio is None else fric_ratio
    friction = tendon.fpj * fr                                      # ΔfpF（摩擦+偏折）
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


# ── 沿長度的預力 P(x) ────────────────────────────────────────────────
# 為什麼要有這一支：compute_losses 是**單點**式——它以一個 fric_ratio 純量代表
# 全長的長度相關損失（摩擦＋滑移），等於「用控制斷面的損失率代表每一斷面」。
# 對 40 m 簡支橋取跨中是合理的（跨中摩擦最大＝最保守，且滑移影響長度未達跨中）；
# 但下列情形會失真，且失真方向**不保守**：
#   ① 錨碇點附近——摩擦≈0 而滑移最大，單點式兩者都取跨中值 → 高估該處 Pe
#   ② 中間錨碇段——P(x) 在錨碇里程處跳升，單點式完全看不到
#   ③ 連續梁——控制斷面不只一個（跨中正彎矩＋墩頂負彎矩），兩處摩擦差很多
#
# loss_profile 於每個里程 x 各自計算，且**用該點的初始力算 f_cgp**：
#     P_i(x) = (f_pj − Δf_fric(x) − Δf_slip(x))·A_ps
# 而 compute_losses 一律用 P_i = f_pj·A_ps（未扣摩擦）。後者高估 f_cgp，
# 因而高估 ES 與潛變——是保守側，但沿長度看會把「錨碇端損失反而較小」這件事抹平。
# 兩式在 μ=K=0、slip=0 時完全相同（見測試 test_loss_profile_degenerates）。


@dataclass
class LossPoint:
    """單一里程的損失分解（應力單位 MPa，Pe 為 N）。"""
    x: float          # 里程 m
    e: float          # 該點偏心 mm（形心下為正）
    fric: float       # 摩擦損失
    slip: float       # 錨具滑移損失
    ES: float         # 彈性縮短
    creep: float      # 潛變
    shrink: float     # 乾縮
    relax: float      # 鬆弛
    short: float      # 短期合計（fric+slip+ES）
    long: float       # 長期合計（creep+shrink+relax）
    total: float      # 總損失
    loss_pct: float   # 總損失率（/fpj）
    fcgp: float       # 鋼腱形心處混凝土應力
    fpe: float        # 有效預力應力
    Pe: float         # 有效預力 N


@dataclass
class LossProfileResult:
    xs: list          # 取樣里程 m
    pts: list         # LossPoint
    Pe_min: float     # 最小有效預力 N
    x_Pemin: float    # 其里程 m
    Pe_max: float
    x_Pemax: float
    Pe_avg: float     # 全長平均（梯形積分／L）
    loss_pct_max: float
    x_lossmax: float
    loss_pct_min: float

    def at(self, x: float) -> LossPoint:
        """任意里程的損失（相鄰取樣點線性內插）。

        ⚠ 內插是線性的，但 f_cgp 沿 x 為二次（M_D 為二次），故取樣點之間有
        內插誤差（n=20、L=40 m 時約 0.02 MPa）。需要精確值時把 n 調大，
        或直接取 pts 中的取樣點。
        """
        xs = self.xs
        if x <= xs[0]:
            return self.pts[0]
        if x >= xs[-1]:
            return self.pts[-1]
        for i in range(len(xs) - 1):
            if xs[i] <= x <= xs[i + 1]:
                t = (x - xs[i]) / (xs[i + 1] - xs[i]) if xs[i + 1] > xs[i] else 0.0
                a, b = self.pts[i], self.pts[i + 1]
                lerp = lambda u, v: u + (v - u) * t
                return LossPoint(
                    x=x, e=lerp(a.e, b.e), fric=lerp(a.fric, b.fric),
                    slip=lerp(a.slip, b.slip), ES=lerp(a.ES, b.ES),
                    creep=lerp(a.creep, b.creep), shrink=a.shrink, relax=a.relax,
                    short=lerp(a.short, b.short), long=lerp(a.long, b.long),
                    total=lerp(a.total, b.total), loss_pct=lerp(a.loss_pct, b.loss_pct),
                    fcgp=lerp(a.fcgp, b.fcgp), fpe=lerp(a.fpe, b.fpe),
                    Pe=lerp(a.Pe, b.Pe))
        return self.pts[-1]


def parabolic_e(e_mid: float, e_end: float, L: float):
    """單跨拋物線偏心 e(x)：端點 e_end、跨中 e_mid，回傳 f(x[m]) → mm。"""
    def f(x: float) -> float:
        return e_end + (e_mid - e_end) * 4.0 * x * (L - x) / (L * L)
    return f


def udl_moment(w_kNm: float, L: float):
    """簡支均佈載重彎矩 M(x)=w·x·(L−x)/2，回傳 f(x[m]) → kN·m。"""
    def f(x: float) -> float:
        return w_kNm * x * (L - x) / 2.0
    return f


def loss_profile(tendon: Tendon, section: Section, L: float,
                 e_fn, MD_fn, segs: list,
                 M_SDL_fn=None, mu: float = 0.25, K: float = 0.003,
                 jack: str = "both", slip_mm: float = 0.0, Ep: float = 195000.0,
                 anchors=None, RH: float = 75.0, relax: float = 10.0,
                 Ep_Eci: float = 7.33, n: int = 20) -> LossProfileResult:
    """沿全長逐點計算預力損失與 P(x)。

    L[m]；e_fn(x[m])→mm（形心下為正）；MD_fn(x)→kN·m（自重，計 ES/潛變用）；
    M_SDL_fn(x)→kN·m（附加恆載，計潛變回復用，不給則為 0）；
    segs：曲率段 [(x1,x2,κ)]（見 tendon_profile.friction_angle）；
    jack/slip_mm/anchors 與 tendon_profile.tendon_forces 同義。

    ⚠ 各點的 f_cgp 以**該點扣除摩擦與滑移後的初始力**計算，與 compute_losses
    以 f_pj·A_ps 計算不同（見本節開頭說明）。
    """
    from .tendon_profile import friction_at, tendon_slip_loss, segmented_tendon_force

    fpj, Aps, N = tendon.fpj, tendon.Aps, tendon.n_tendons
    shrink = 0.8 * (1195 - 10.55 * RH) * 0.0981
    xs = [L * i / n for i in range(n + 1)]
    pts = []
    for x in xs:
        if anchors:
            f = segmented_tendon_force(x, anchors, segs, fpj, Aps, 0.0,
                                       mu, K, jack, slip_mm, Ep)
            fric, sl = f.friction, f.slip
        else:
            fric = fpj * friction_at(x, L, segs, mu, K, jack)
            sl = tendon_slip_loss(x, L, segs, fpj, jack, mu, K, slip_mm, Ep)
        Pi_x = (fpj - fric - sl) * Aps
        e = e_fn(x)
        M_D = MD_fn(x) * 1e6
        M_SDL = (M_SDL_fn(x) * 1e6) if M_SDL_fn else 0.0
        fcgp = Pi_x / section.A + Pi_x * e * e / section.I - M_D * e / section.I
        fcds = M_SDL * e / section.I
        ES = (N - 1) / (2 * N) * Ep_Eci * fcgp
        creep = 12.0 * fcgp - 7.0 * fcds
        short = fric + sl + ES
        lng = creep + shrink + relax
        total = short + lng
        fpe = fpj - total
        pts.append(LossPoint(x=x, e=e, fric=fric, slip=sl, ES=ES, creep=creep,
                             shrink=shrink, relax=relax, short=short, long=lng,
                             total=total, loss_pct=total / fpj, fcgp=fcgp,
                             fpe=fpe, Pe=fpe * Aps))
    Pes = [p.Pe for p in pts]
    lps = [p.loss_pct for p in pts]
    imin, imax = Pes.index(min(Pes)), Pes.index(max(Pes))
    ilmax = lps.index(max(lps))
    area = sum((Pes[i] + Pes[i + 1]) / 2 * (xs[i + 1] - xs[i]) for i in range(n))
    return LossProfileResult(
        xs=xs, pts=pts, Pe_min=Pes[imin], x_Pemin=xs[imin],
        Pe_max=Pes[imax], x_Pemax=xs[imax], Pe_avg=area / L if L else 0.0,
        loss_pct_max=lps[ilmax], x_lossmax=xs[ilmax], loss_pct_min=min(lps))
