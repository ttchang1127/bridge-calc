"""連續梁變斷面（中墩底板加厚）：斷面沿長度變化、EI 變化的柔度與支承彎矩。

對應 算例_連續梁次彎矩與服務性應力驗算、analyzer ⑦ 中墩形心回饋。
單位：長度 m（沿橋）、mm（斷面）、載重 kN、彎矩 kN·m（正彎矩為正）。

作法（力法，EI(x)＝EI_ref·I_rel(x)）：
  柔度 F_ij ＝ ∫ m_i·m_j / I_rel dx
  單位集中載重於 p 的轉角項 δ_i(p) ＝ ∫ M0(x;p)·m_i(x)/I_rel dx
      依 Maxwell 互易＝「簡支跨以曲率 κ_i = m_i/I_rel 產生的撓度」u_i 在 p 的值：
      u_i″ = −κ_i、u_i(0)=u_i(L)=0 → u_i(t) = −S2(t) + t·S2(L)/L（S1=∫κ、S2=∫S1）
      故每個支承只需一次累積積分，之後任意載重位置 O(1) 查表（三次 Hermite 內插）。
  支承彎矩 X ＝ −F⁻¹·δ。等斷面時 δ_i(p) 退化為 P·a·b·(L+a)/(6L)（與三彎矩閉合解同）。

⚠ 偏心參考軸仍為等斷面（基準）形心：鋼腱絕對高度 y = ȳb_ref − e_ref；
   相對局部形心的偏心 e_loc(x) = e_ref − Δȳ(x)，Δȳ = ȳb_ref − ȳb(x)（加厚→形心下移→Δȳ>0）。
"""
from dataclasses import dataclass, field
from math import sqrt
from typing import Callable, List, Optional, Sequence

from .model import Section

_GL = ((-sqrt(3 / 5), 5 / 9), (0.0, 8 / 9), (sqrt(3 / 5), 5 / 9))   # 3 點 Gauss-Legendre（五次精確）


def section_from_dims(top_w: float, top_t: float, bot_w: float, bot_t: float,
                      web_t: float, n_web: int, h: float) -> Section:
    """單室箱梁由構件尺寸算斷面性質（頂板/底板/腹板三塊矩形；同 engine.js BC.sectionFromDims）。"""
    hw = h - top_t - bot_t
    parts = [(top_w * top_t, h - top_t / 2, top_w * top_t ** 3 / 12),
             (bot_w * bot_t, bot_t / 2, bot_w * bot_t ** 3 / 12),
             (n_web * web_t * hw, bot_t + hw / 2, n_web * web_t * hw ** 3 / 12)]
    A = sum(p[0] for p in parts)
    yb = sum(p[0] * p[1] for p in parts) / A
    I = sum(p[2] + p[0] * (p[1] - yb) ** 2 for p in parts)
    return Section(A, I, yb, h)


@dataclass
class HaunchProfile:
    """中墩底板加厚：距各內支承 d < length 時底板厚由 bot_t 線性增至 bot_t_pier。"""
    spans: List[float]
    xs: List[float]
    dims: tuple                # (top_w, top_t, bot_w, bot_t, web_t, n_web, h)
    bot_t_pier: float
    lengths: List[float]       # 各內支承採用的單側加厚長度（≤0.5×相鄰跨）
    sec_ref: Section
    breaks: List[float]
    _cache: dict = field(default_factory=dict, repr=False)

    def bot_t_at(self, x: float) -> float:
        top_w, top_t, bot_w, bot_t, web_t, n_web, h = self.dims
        t = bot_t
        for j in range(1, len(self.xs) - 1):
            d = abs(x - self.xs[j])
            b = self.lengths[j - 1]
            if b > 0 and d < b:
                t = max(t, self.bot_t_pier + (bot_t - self.bot_t_pier) * d / b)
        return t

    def section_at(self, x: float) -> Section:
        key = round(x, 6)
        s = self._cache.get(key)
        if s is None:
            top_w, top_t, bot_w, bot_t, web_t, n_web, h = self.dims
            s = self._cache[key] = section_from_dims(top_w, top_t, bot_w, self.bot_t_at(x), web_t, n_web, h)
        return s

    def I_rel(self, x: float) -> float:
        return self.section_at(x).I / self.sec_ref.I

    def A_rel(self, x: float) -> float:
        return self.section_at(x).A / self.sec_ref.A

    def dyb(self, x: float) -> float:
        """形心下移量 Δȳ = ȳb_ref − ȳb(x)（mm）。"""
        return self.sec_ref.yb - self.section_at(x).yb


def haunch_profile(spans: Sequence[float], top_w: float, top_t: float, bot_w: float, bot_t: float,
                   web_t: float, n_web: int, h: float, bot_t_pier: float, length: float) -> HaunchProfile:
    xs = [0.0]
    for L in spans:
        xs.append(xs[-1] + L)
    lens, brk = [], []
    for j in range(1, len(spans)):
        b = max(0.0, min(length, 0.5 * spans[j - 1], 0.5 * spans[j]))
        lens.append(b)
        brk += [xs[j] - b, xs[j], xs[j] + b]
    dims = (top_w, top_t, bot_w, bot_t, web_t, n_web, h)
    return HaunchProfile(list(spans), xs, dims, max(bot_t, bot_t_pier), lens,
                         section_from_dims(*dims), sorted(set(round(v, 9) for v in brk)))


def _cells(a: float, b: float, cuts: Sequence[float], h: float) -> List[tuple]:
    pts = sorted(set([a, b] + [c for c in cuts if a < c < b]))
    out = []
    for u, v in zip(pts[:-1], pts[1:]):
        n = max(1, int(round((v - u) / h)))
        d = (v - u) / n
        out += [(u + k * d, u + (k + 1) * d) for k in range(n)]
    return out


def gauss_integrate(f: Callable[[float], float], a: float, b: float,
                    cuts: Sequence[float] = (), h: float = 0.25) -> float:
    """∫f，於 cuts 處切開（容許跳躍）、每小段 3 點 Gauss。"""
    s = 0.0
    for u, v in _cells(a, b, cuts, h):
        m, r = 0.5 * (u + v), 0.5 * (v - u)
        s += r * sum(w * f(m + r * g) for g, w in _GL)
    return s


class ContFlex:
    """連續梁柔度（變 EI）：F、各內支承的 u_i 撓度表，供點載重／分佈載重求支承彎矩。"""

    def __init__(self, spans: Sequence[float], I_rel: Callable[[float], float],
                 breaks: Sequence[float] = (), grid: float = 0.05):
        self.spans = list(spans)
        xs = [0.0]
        for L in spans:
            xs.append(xs[-1] + L)
        self.xs, self.I_rel, self.breaks = xs, I_rel, list(breaks)
        n = len(spans)
        self.nu = nu = n - 1
        self.F = [[0.0] * nu for _ in range(nu)]
        cuts = xs + self.breaks
        for i in range(nu):                      # 內支承 i+1 的單位三角形 m 非零於 [xs[i], xs[i+2]]
            self.F[i][i] = gauss_integrate(lambda x: self._m(i + 1, x) ** 2 / I_rel(x), xs[i], xs[i + 2], cuts, 0.25)
            if i + 1 < nu:                       # 與內支承 i+2 重疊於跨 [xs[i+1], xs[i+2]]
                self.F[i][i + 1] = self.F[i + 1][i] = gauss_integrate(
                    lambda x: self._m(i + 1, x) * self._m(i + 2, x) / I_rel(x), xs[i + 1], xs[i + 2], cuts, 0.25)
        # 各內支承 i：左跨（i−1）與右跨（i）的撓度表 u、u′
        self.tab = {}
        for i in range(1, n):
            for k in (i - 1, i):
                self.tab[(i, k)] = self._deflection_table(i, k, grid)

    def _m(self, i: int, x: float) -> float:
        xs = self.xs
        a, c, b = xs[i - 1], xs[i], xs[i + 1]
        if a - 1e-12 <= x <= c:
            return (x - a) / (c - a)
        if c <= x <= b + 1e-12:
            return (b - x) / (b - c)
        return 0.0

    def _deflection_table(self, i: int, k: int, grid: float):
        x0, L = self.xs[k], self.spans[k]
        cells = _cells(x0, x0 + L, self.breaks, grid)
        ts, S1, S2 = [x0], [0.0], [0.0]
        for u, v in cells:
            m, r = 0.5 * (u + v), 0.5 * (v - u)
            k1 = k2 = 0.0
            for g, w in _GL:
                t = m + r * g
                kap = self._m(i, t) / self.I_rel(t)
                k1 += r * w * kap
                k2 += r * w * kap * (v - t)
            S2.append(S2[-1] + (v - u) * S1[-1] + k2)
            S1.append(S1[-1] + k1)
            ts.append(v)
        c = S2[-1] / L
        u_ = [-s2 + (t - x0) * c for t, s2 in zip(ts, S2)]
        du = [-s1 + c for s1 in S1]
        return ts, u_, du

    def _u(self, i: int, k: int, p: float) -> float:
        ts, u, du = self.tab[(i, k)]
        lo, hi = 0, len(ts) - 1
        if p <= ts[0]:
            return u[0]
        if p >= ts[-1]:
            return u[-1]
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if ts[mid] <= p:
                lo = mid
            else:
                hi = mid
        t0, t1 = ts[lo], ts[hi]
        hh = t1 - t0
        s = (p - t0) / hh
        h00, h10 = 2 * s ** 3 - 3 * s ** 2 + 1, s ** 3 - 2 * s ** 2 + s
        h01, h11 = -2 * s ** 3 + 3 * s ** 2, s ** 3 - s ** 2
        return h00 * u[lo] + h10 * hh * du[lo] + h01 * u[hi] + h11 * hh * du[hi]

    def _solve(self, delta: List[float]) -> List[float]:
        nu = self.nu
        if nu <= 0:
            return [0.0] * (len(self.spans) + 1)
        A = [row[:] + [-delta[i]] for i, row in enumerate(self.F)]
        for i in range(nu):
            for j in range(i + 1, nu):
                f = A[j][i] / A[i][i]
                for c in range(i, nu + 1):
                    A[j][c] -= f * A[i][c]
        X = [0.0] * nu
        for i in range(nu - 1, -1, -1):
            X[i] = (A[i][nu] - sum(A[i][c] * X[c] for c in range(i + 1, nu))) / A[i][i]
        return [0.0] + X + [0.0]

    def span_of(self, x: float) -> int:
        for i in range(len(self.spans)):
            if x <= self.xs[i + 1] + 1e-9:
                return i
        return len(self.spans) - 1

    def support_moments_point(self, p: float, P: float = 1.0) -> List[float]:
        n = len(self.spans)
        if n < 2 or p < -1e-9 or p > self.xs[-1] + 1e-9:
            return [0.0] * (n + 1)
        k = self.span_of(p)
        delta = [0.0] * self.nu
        for i in (k, k + 1):
            if 1 <= i <= n - 1:
                delta[i - 1] = P * self._u(i, k, p)
        return self._solve(delta)

    def support_moments_dist(self, w: Callable[[float], float]) -> List[float]:
        n = len(self.spans)
        delta = [0.0] * self.nu
        for i in range(1, n):
            for k in (i - 1, i):
                delta[i - 1] += gauss_integrate(lambda t: w(t) * self._u(i, k, t),
                                                self.xs[k], self.xs[k + 1], self.breaks, 0.25)
        return self._solve(delta)

    def moment_at(self, X: List[float], x: float, M0: float) -> float:
        i = self.span_of(x)
        t = (x - self.xs[i]) / self.spans[i]
        return M0 + X[i] * (1 - t) + X[i + 1] * t

    def M0_dist(self, w: Callable[[float], float], x: float) -> float:
        """簡支跨以分佈載重 w(t) 於斷面 x 的彎矩 ∫G(x,t)w(t)dt。"""
        k = self.span_of(x)
        x0, L = self.xs[k], self.spans[k]
        a = x - x0

        def g(t):
            s = t - x0
            return (s * (L - a) / L if s <= a else a * (L - s) / L) * w(t)
        return gauss_integrate(g, x0, x0 + L, self.breaks + [x], 0.25)

    def dl_moment(self, w: Callable[[float], float], x: float, X: Optional[List[float]] = None) -> float:
        X = X if X is not None else self.support_moments_dist(w)
        return self.moment_at(X, x, self.M0_dist(w, x))
