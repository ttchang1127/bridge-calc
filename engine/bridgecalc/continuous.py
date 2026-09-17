"""連續梁中墩（次彎矩 M2、T 斷面極限彎曲、墩斷面服務性）。

對應 算例_連續梁次彎矩與服務性應力驗算、算例_極限強度彎曲設計 §十二（中墩）。
40+40 兩跨連續後張箱梁（2026-09-17 校正：M2 改力法實算、頂板腱 e=646）：
  - M2 全長正彎矩：跨中不利、墩頂有利（原算例正負號相反）
  - T 斷面 M1：負彎矩使 NA 進腹板（c>hf）→ Mn 縮水 → CR≈0.55 嚴重不足（Mu 已計入 M2）
  - B 墩服務性：含 M2 底緣 −12.70 MPa ✓；不計 M2 為 −19.72（假性超限）
單位：力 N、長度 mm、應力 MPa；M 以 kN·m 表示。
"""
from dataclasses import dataclass
from math import sqrt
from typing import Callable, List, Optional, Sequence, Union
from .model import Section
from .flexure import beta1


def secondary_moment(M_total_kNm: float, M1_kNm: float) -> float:
    """次彎矩 M2 = M_total − M1（超靜定連續梁；M1=ΣPe·e 一次彎矩）。"""
    return M_total_kNm - M1_kNm


def primary_moment(layers) -> float:
    """一次預力彎矩 M1 = Σ(Pe_layer · e_layer)。layers: [(Pe_kN, e_m), ...]。"""
    return sum(Pe * e for Pe, e in layers)


@dataclass
class TFlexureResult:
    c: float          # 中性軸 mm
    flanged: bool     # True=T 斷面（NA 進腹板）
    a: float          # 等效矩形深 mm
    fps: float        # 鋼腱應力 MPa
    Mn: float         # 公稱彎矩 kN·m
    eps_t: float      # 極限拉應變
    phi: float        # 強度折減
    phiMn: float      # 設計彎矩 kN·m
    CR: float         # 強度比 φMn/Mu
    ok: bool


def flexural_strength_T(Aps: float, fpu: float, fc: float,
                        b: float, hf: float, bw: float, dp: float,
                        Mu_kNm: float, dt: float = None) -> TFlexureResult:
    """T／矩形斷面極限彎曲：先試矩形，c>hf 改 T 斷面（壓力區進腹板）。

    b 受壓翼緣有效寬、hf 翼緣厚、bw 腹板寬、dp 鋼腱有效深。
    中墩負彎矩 → 受壓區在底板（窄）→ 常進腹板需 T 斷面公式。
    """
    fpy = 0.90 * fpu
    k = 2 * (1.04 - fpy / fpu)
    b1 = beta1(fc)
    c_rect = Aps * fpu / (0.85 * fc * b1 * b + k * Aps * fpu / dp)
    if c_rect <= hf:
        flanged = False
        c = c_rect
        a = b1 * c
        fps = fpu * (1 - k * c / dp)
        Mn = Aps * fps * (dp - a / 2) / 1e6
    else:
        flanged = True
        c = (Aps * fpu - 0.85 * fc * (b - bw) * hf) / (0.85 * fc * b1 * bw + k * Aps * fpu / dp)
        a = b1 * c
        fps = fpu * (1 - k * c / dp)
        Mn = (Aps * fps * (dp - a / 2)
              + 0.85 * fc * (b - bw) * hf * (a / 2 - hf / 2)) / 1e6
    dt = dt if dt is not None else dp
    eps_t = (dt - c) / c * 0.003
    if eps_t >= 0.005:
        phi = 1.0
    elif eps_t <= 0.002:
        phi = 0.75
    else:
        phi = 0.75 + 0.25 * (eps_t - 0.002) / 0.003
    phiMn = phi * Mn
    return TFlexureResult(c, flanged, a, fps, Mn, eps_t, phi, phiMn,
                          phiMn / Mu_kNm, phiMn >= Mu_kNm)


def pier_service_stress(Pe: float, section: Section, e: float, M_ext_kNm: float):
    """連續梁墩斷面服務性應力（壓為負）。

    e 慣例同 service.stresses()：**形心下為正**；負彎矩區頂板 PT 在形心上 → e 取負。
    M_ext 為外力彎矩（含 M2 = M_DL+M_LL+M2），負彎矩（hogging）取負值。
    σ_t = −Pe/A + Pe·e/St − M_ext/St ；σ_b = −Pe/A − Pe·e/Sb + M_ext/Sb。回傳 (σ_t, σ_b)。
    例：B 墩 Pe=36,257kN、e=−276（頂板PT合力形心上）、M_ext=−41,080+M2 17,354 → σ_b=−12.70 MPa。
    """
    M = M_ext_kNm * 1e6
    st = -Pe / section.A + Pe * e / section.St - M / section.St
    sb = -Pe / section.A - Pe * e / section.Sb + M / section.Sb
    return st, sb


# ════════════════ 連續梁鋼腱線形（分段拋物線）＋ 次彎矩 M2（力法）════════════════
# 慣例：x 以 m、e 以 mm（**形心下為正**）、P 以 kN、M 以 kN·m（**正彎矩＝底緣受拉為正**）。
# 一次彎矩 M1 = −P·e/1000（鋼腱在形心下 → 上拱 → 負彎矩）。
@dataclass
class ParabolaSeg:
    """頂點式拋物線段 e(x) = ev + c·(x − xv)²，適用 x1 ≤ x ≤ x2。c 單位 mm/m²。"""
    x1: float
    x2: float
    xv: float
    ev: float
    c: float
    kind: str = ""
    R: float = float("inf")   # 曲率半徑 mm = 10⁶/(2|c|)

    def e(self, x: float) -> float:
        d = x - self.xv
        return self.ev + self.c * d * d

    def slope(self, x: float) -> float:
        """de/dx（mm/m）。"""
        return 2 * self.c * (x - self.xv)


def parabola_seg(x1, x2, xv, ev, c, kind="") -> ParabolaSeg:
    R = 1e6 / (2 * abs(c)) if abs(c) > 1e-12 else float("inf")
    return ParabolaSeg(x1, x2, xv, ev, c, kind, R)


@dataclass
class ContTendonProfile:
    segs: List[ParabolaSeg]
    sx: List[float]        # 支承里程 m
    lows: List[float]      # 各跨低點里程 m
    infl: List[float]      # 反曲點里程 m
    total: float
    k_pier: float
    r_end: float
    R_min: float           # 各段最小曲率半徑 mm

    def e_at(self, x: float) -> float:
        for g in self.segs:
            if g.x1 - 1e-9 <= x <= g.x2 + 1e-9:
                return g.e(x)
        return self.segs[-1].e(x)


def cont_tendon_segs(spans: Sequence[float], e_end: float, e_mid: float, e_pier: float,
                     k_pier: float = 0.12, r_end: float = 0.42) -> ContTendonProfile:
    """連續梁全長連續腱線形（G1 §五之補、FHWA Eqn. 3.30–3.35）。

    墩頂段頂點在墩心（斜率 0，長 b₁ = k_pier×該側跨長）；跨中段頂點在低點
    （端跨 r_end×跨長、內跨中點）；兩段坡度連續 h₁/b₁ = h₂/b₂ → h_i = Δ·b_i/s，
    反曲點落在交界。端跨端錨→低點為單段、頂點在低點。
    """
    n = len(spans)
    xs = [0.0]
    for L in spans:
        xs.append(xs[-1] + L)
    lows = []
    for i, L in enumerate(spans):
        if n == 1:
            lows.append(xs[0] + L / 2)
        elif i == 0:
            lows.append(xs[0] + r_end * L)
        elif i == n - 1:
            lows.append(xs[n - 1] + (1 - r_end) * L)
        else:
            lows.append(xs[i] + L / 2)
    segs: List[ParabolaSeg] = []
    infl: List[float] = []

    def push(x1, x2, xv, ev, c, kind):
        if x2 > x1 + 1e-9:
            segs.append(parabola_seg(x1, x2, xv, ev, c, kind))

    def half_to_pier(xl, xp, L, left_to_right):
        s = abs(xp - xl)
        b1 = min(k_pier * L, s * 0.9)
        b2 = s - b1
        D = e_mid - e_pier
        h1, h2 = D * b1 / s, D * b2 / s
        if left_to_right:
            push(xl, xp - b1, xl, e_mid, -h2 / (b2 * b2), "sag")
            push(xp - b1, xp, xp, e_pier, h1 / (b1 * b1), "hog")
            infl.append(xp - b1)
        else:
            push(xp, xp + b1, xp, e_pier, h1 / (b1 * b1), "hog")
            push(xp + b1, xl, xl, e_mid, -h2 / (b2 * b2), "sag")
            infl.append(xp + b1)

    for i in range(n):
        L, xa, xb, xl = spans[i], xs[i], xs[i + 1], lows[i]
        if i == 0:
            be, he = xl - xa, e_mid - e_end
            push(xa, xl, xl, e_mid, -he / (be * be), "sag")
        else:
            half_to_pier(xl, xa, L, False)
        if i == n - 1:
            be, he = xb - xl, e_mid - e_end
            push(xl, xb, xl, e_mid, -he / (be * be), "sag")
        else:
            half_to_pier(xl, xb, L, True)
    return ContTendonProfile(segs, xs, lows, infl, xs[-1], k_pier, r_end,
                             min(g.R for g in segs))


@dataclass
class TendonGroup:
    """一組鋼腱（同線形）：P 為常數 kN 或 P(x) 函式；segs 未覆蓋處＝該組不存在（如僅配墩頂的頂板腱）。"""
    P: Union[float, Callable[[float], float]]
    segs: List[ParabolaSeg]


def _group_seg(g: TendonGroup, x: float, probe: Optional[float] = None) -> Optional[ParabolaSeg]:
    t = x if probe is None else probe
    for s in g.segs:
        if s.x1 - 1e-9 <= t <= s.x2 + 1e-9:
            return s
    return None


def primary_moment_at(groups: Sequence[TendonGroup], x: float, probe: Optional[float] = None) -> float:
    """一次彎矩 M1(x) = Σ −P·e/1000（kN·m，正彎矩為正）。

    probe：用來判定「哪一段／該組是否存在」的里程（積分時取子區間中點，
    使錨碇處 M1 的跳躍取到正確的單側極限）。
    """
    M = 0.0
    for g in groups:
        s = _group_seg(g, x, probe)
        if s is None:
            continue
        P = g.P(x) if callable(g.P) else g.P
        M += -P * s.e(x) / 1000.0
    return M


def group_breaks(groups: Sequence[TendonGroup]) -> List[float]:
    out = set()
    for g in groups:
        for s in g.segs:
            out.add(round(s.x1, 9))
            out.add(round(s.x2, 9))
    return sorted(out)


@dataclass
class ForceMethodM2Result:
    sx: List[float]
    X: List[float]         # 各支承 M2（kN·m），兩端恆為 0
    F: List[List[float]]   # 柔度係數 ∫m_i m_j dx（EI 常數）
    b: List[float]         # ∫M1·m_i dx

    def M2_at(self, x: float) -> float:
        sx = self.sx
        for i in range(len(sx) - 1):
            if sx[i] - 1e-9 <= x <= sx[i + 1] + 1e-9:
                t = (x - sx[i]) / (sx[i + 1] - sx[i])
                return self.X[i] * (1 - t) + self.X[i + 1] * t
        return 0.0


def secondary_moments_force(spans: Sequence[float],
                            M1: Callable[[float, float], float],
                            breaks: Sequence[float] = (), n_sub: int = 8) -> ForceMethodM2Result:
    """連續梁預力次彎矩 M2（力法／柔度法，EI 常數）。

    以各內支承彎矩 X_i 為贅餘力，M2 在支承間為直線、端支承為 0：
        M2(x) = Σ X_i·m_i(x)，m_i ＝支承 i 處為 1、相鄰支承為 0 的三角形
    轉角連續（相容）條件：∫(M1 + M2)·m_i dx = 0 → Σ_j F_ij X_j = −∫M1·m_i dx，
        F_ii = (L_左 + L_右)/3、F_i,i+1 = L/6（即三彎矩方程，M1 當「曲率載重」）。

    M1(x, probe)：probe 為所在子區間中點，供錨碇跳躍取單側極限（見 primary_moment_at）。
    breaks：M1 的折點／跳躍里程（線形段交界、錨碇點）；支承里程自動加入。
    積分以 Simpson，子區間對齊所有折點；M1 為分段二次、m_i 分段線性 → 被積函數分段三次，
    Simpson **精確**（P 沿長度變時 n_sub 加密）。

    與等效載重法的關係：兩者等價，但力法**不需**把 M1 二次微分成等效載重，
    故端錨偏心彎矩、墩頂折角集中力、部分長度鋼腱的錨碇力都自動計入。
    """
    n = len(spans)
    sx = [0.0]
    for L in spans:
        sx.append(sx[-1] + L)
    if n < 2:
        return ForceMethodM2Result(sx, [0.0] * (n + 1), [], [])
    pts = sorted(set([round(v, 9) for v in sx] + [round(v, 9) for v in breaks
                                                  if sx[0] < v < sx[-1]]))

    def m_hat(i, x):          # 內支承 i（1..n−1）的單位三角形
        a, c0, b = sx[i - 1], sx[i], sx[i + 1]
        if a <= x <= c0:
            return (x - a) / (c0 - a)
        if c0 <= x <= b:
            return (b - x) / (b - c0)
        return 0.0

    nu = n - 1
    bvec = [0.0] * nu
    for k in range(len(pts) - 1):
        a, b = pts[k], pts[k + 1]
        if b - a < 1e-9:
            continue
        mid = 0.5 * (a + b)
        h = (b - a) / (2 * n_sub)
        for i in range(nu):
            if not (sx[i] - 1e-9 <= mid <= sx[i + 2] + 1e-9):
                continue
            s = 0.0
            for j in range(2 * n_sub + 1):
                x = a + j * h
                w = 1 if j in (0, 2 * n_sub) else (4 if j % 2 else 2)
                s += w * M1(x, mid) * m_hat(i + 1, x)
            bvec[i] += s * h / 3
    F = [[0.0] * nu for _ in range(nu)]
    for i in range(nu):
        F[i][i] = (spans[i] + spans[i + 1]) / 3
        if i + 1 < nu:
            F[i][i + 1] = F[i + 1][i] = spans[i + 1] / 6
    # 高斯消去（F 為對稱正定三對角）
    A = [row[:] + [-bvec[i]] for i, row in enumerate(F)]
    for i in range(nu):
        p = A[i][i]
        for j in range(i + 1, nu):
            f = A[j][i] / p
            for k in range(i, nu + 1):
                A[j][k] -= f * A[i][k]
    Xi = [0.0] * nu
    for i in range(nu - 1, -1, -1):
        Xi[i] = (A[i][nu] - sum(A[i][k] * Xi[k] for k in range(i + 1, nu))) / A[i][i]
    return ForceMethodM2Result(sx, [0.0] + Xi + [0.0], F, bvec)


def continuous_prestress(spans: Sequence[float], groups: Sequence[TendonGroup],
                         n_sub: int = 8) -> ForceMethodM2Result:
    """便利函式：由鋼腱組直接求 M2（M1 = Σ −P·e/1000）。"""
    return secondary_moments_force(
        spans, lambda x, probe: primary_moment_at(groups, x, probe),
        group_breaks(groups), n_sub)
