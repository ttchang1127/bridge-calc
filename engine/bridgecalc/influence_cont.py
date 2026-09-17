"""連續梁影響線與台灣 HS20-44 設計彎矩包絡（解析三彎矩法，EI 常數）。

對應 活載計算詳程快查卡 §二（連續梁佈置）/§四（衝擊長度）、算例_連續梁次彎矩與服務性應力驗算 §六。
單位：長度 m、載重 kN → 彎矩 kN·m（正彎矩為正）。

作法：單位集中載重於 p 時，內支承彎矩 X 由三彎矩方程（力法）解出：
    L_i·X_{i−1} + 2(L_i+L_{i+1})·X_i + L_{i+1}·X_{i+1} = −6·∫M0·m_i dx
    點載重（距該跨左端 a、b=L−a）：∫M0·m = P·a·b·(L+a)/(6L)（對右端支承）、P·a·b·(L+b)/(6L)（對左端）
斷面 x 的彎矩 = 簡支彎矩 M0（載重與斷面同跨時）＋ 支承彎矩在該跨的線性內插。
影響線值為解析解（無有限元素節點內插誤差）；車道均布面積以細格梯形分正負積分（跨零點切開）。

台灣規範（§3.8、§3.9、§3.13）：
  - 設計卡車 HS20-44 與車道載重「取大」，每車道一輛車；車可雙向行駛
  - 車道載重：均布 9.4 kN/m 載於同號影響線區＋集中載重 80 kN（彎矩）
  - **最大負彎矩**：另加一個相等集中載重於其他跨（共 2 個，放於不同跨之最負縱距）
  - 衝擊 I = 15.24/(L+38.1) ≤ 0.30；L：正彎矩取該跨跨徑、**負彎矩取相鄰兩載重跨平均**
  - 多車道折減：1~2 車道 1.00、3 車道 0.90、≥4 車道 0.75
"""
from dataclasses import dataclass
from typing import List, Sequence

from .influence import (TW_HS20_AXLES, TW_HS20_SPACING, TW_LANE_W, TW_LANE_PM,
                        taiwan_impact)


def _supports(spans: Sequence[float]) -> List[float]:
    xs = [0.0]
    for L in spans:
        xs.append(xs[-1] + L)
    return xs


def _span_of(xs: Sequence[float], x: float) -> int:
    n = len(xs) - 1
    for i in range(n):
        if x <= xs[i + 1] + 1e-9:
            return i
    return n - 1


def _solve_tridiag(spans: Sequence[float], rhs: Sequence[float]) -> List[float]:
    """解三彎矩方程（內支承 1..n−1），回傳含兩端 0 的支承彎矩。rhs[i] 對應內支承 i+1。"""
    n = len(spans)
    nu = n - 1
    if nu <= 0:
        return [0.0] * (n + 1)
    a = [spans[i] for i in range(nu)]            # 下對角（X_{i−1} 係數，i≥1 才用）
    b = [2 * (spans[i] + spans[i + 1]) for i in range(nu)]
    c = [spans[i + 1] for i in range(nu)]        # 上對角
    d = list(rhs)
    for i in range(1, nu):                       # Thomas
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    X = [0.0] * nu
    X[-1] = d[-1] / b[-1]
    for i in range(nu - 2, -1, -1):
        X[i] = (d[i] - c[i] * X[i + 1]) / b[i]
    return [0.0] + X + [0.0]


def cont_support_moments_point(spans: Sequence[float], p: float, P: float = 1.0) -> List[float]:
    """集中載重 P 於里程 p 時各支承彎矩（kN·m，兩端 0）。"""
    xs = _supports(spans)
    n = len(spans)
    rhs = [0.0] * (n - 1)
    if n < 2 or p < -1e-9 or p > xs[-1] + 1e-9:
        return [0.0] * (n + 1)
    k = _span_of(xs, p)
    L = spans[k]
    aa = min(max(p - xs[k], 0.0), L)
    bb = L - aa
    if k >= 1:                                   # 本跨左端＝內支承 k
        rhs[k - 1] += -6 * P * aa * bb * (L + bb) / (6 * L)
    if k <= n - 2:                               # 本跨右端＝內支承 k+1
        rhs[k] += -6 * P * aa * bb * (L + aa) / (6 * L)
    return _solve_tridiag(spans, rhs)


def cont_support_moments_uniform(spans: Sequence[float], w: Sequence[float]) -> List[float]:
    """各跨均布載重 w[i]（kN/m）時各支承彎矩；∫M0·m = w·L³/24（兩端相同）。"""
    n = len(spans)
    rhs = [0.0] * (n - 1)
    for k in range(n):
        t = -6 * w[k] * spans[k] ** 3 / 24
        if k >= 1:
            rhs[k - 1] += t
        if k <= n - 2:
            rhs[k] += t
    return _solve_tridiag(spans, rhs)


def _moment_at(spans, xs, X, x, M0) -> float:
    i = _span_of(xs, x)
    t = (x - xs[i]) / spans[i]
    return M0 + X[i] * (1 - t) + X[i + 1] * t


def cont_moment_il(spans: Sequence[float], x: float, p: float) -> float:
    """斷面 x 的彎矩影響線縱距（單位載重於 p）。"""
    xs = _supports(spans)
    X = cont_support_moments_point(spans, p)
    i, k = _span_of(xs, x), _span_of(xs, p)
    M0 = 0.0
    if i == k:
        L = spans[i]
        a, q = x - xs[i], p - xs[i]
        M0 = q * (L - a) / L if q <= a else a * (L - q) / L
    return _moment_at(spans, xs, X, x, M0)


def cont_dl_moment(spans: Sequence[float], w: float, x: float) -> float:
    """全長均布載重 w 於斷面 x 的彎矩（解析）。"""
    xs = _supports(spans)
    X = cont_support_moments_uniform(spans, [w] * len(spans))
    i = _span_of(xs, x)
    a = x - xs[i]
    M0 = w * a * (spans[i] - a) / 2
    return _moment_at(spans, xs, X, x, M0)


def taiwan_lane_reduction(lanes: int) -> float:
    """台灣多車道折減：1~2 車道 1.00、3 車道 0.90、≥4 車道 0.75。"""
    return 1.0 if lanes <= 2 else (0.90 if lanes == 3 else 0.75)


def cont_impact_length(spans: Sequence[float], x: float, sign: int) -> float:
    """衝擊長度：正彎矩＝該跨跨徑；負彎矩＝相鄰兩載重跨平均（取斷面最近內支承的兩側跨）。"""
    xs = _supports(spans)
    n = len(spans)
    i = _span_of(xs, x)
    if sign > 0 or n == 1:
        return spans[i]
    for j in range(1, n):                        # 恰在內支承
        if abs(x - xs[j]) < 1e-6:
            return (spans[j - 1] + spans[j]) / 2
    left_int = i >= 1
    right_int = i <= n - 2
    if left_int and (not right_int or x - xs[i] <= xs[i + 1] - x):
        return (spans[i - 1] + spans[i]) / 2
    return (spans[i] + spans[i + 1]) / 2


@dataclass
class ContLiveMoment:
    x: float
    truck_pos: float
    truck_neg: float
    lane_pos: float
    lane_neg: float
    pos: float          # 每車道、不含衝擊，取大
    neg: float
    I_pos: float
    I_neg: float


class _ILCache:
    def __init__(self, spans):
        self.spans = list(spans)
        self.xs = _supports(spans)
        self.tot = self.xs[-1]
        self._X = {}

    def X(self, p):
        key = round(p, 6)
        v = self._X.get(key)
        if v is None:
            v = self._X[key] = cont_support_moments_point(self.spans, p)
        return v

    def eta(self, x, p):
        xs, spans = self.xs, self.spans
        X = self.X(p)
        i, k = _span_of(xs, x), _span_of(xs, p)
        M0 = 0.0
        if i == k:
            L = spans[i]
            a, q = x - xs[i], p - xs[i]
            M0 = q * (L - a) / L if q <= a else a * (L - q) / L
        return _moment_at(spans, xs, X, x, M0)


def _live_at(c: _ILCache, x: float, step: float, grid_per_span: int) -> ContLiveMoment:
    spans, xs, tot = c.spans, c.xs, c.tot
    # 設計卡車：雙向、整數步進、軸位夾到 [0, tot]
    span_t = TW_HS20_SPACING[-1]
    sets = ((TW_HS20_AXLES, TW_HS20_SPACING),
            (tuple(reversed(TW_HS20_AXLES)), tuple(span_t - d for d in reversed(TW_HS20_SPACING))))
    nstep = int(round((tot + span_t) / step))
    tp = tn = 0.0
    for P_set, d_set in sets:
        for k in range(nstep + 1):
            s = -span_t + k * step
            v = 0.0
            for P, d in zip(P_set, d_set):
                ax = s + d
                if -1e-9 <= ax <= tot + 1e-9:
                    v += P * c.eta(x, min(max(ax, 0.0), tot))
            tp, tn = max(tp, v), min(tn, v)
    # 車道：分跨細格、逐小段依正負拆分（跨零點以線性內插切開）積分同號面積；
    # 集中載重取最大縱距，負彎矩另加他跨一個（共 2 個，取不同跨之最負縱距）
    pos_area = neg_area = 0.0
    span_min = []
    eta_max = 0.0
    for j, L in enumerate(spans):
        h = L / grid_per_span
        vals = [c.eta(x, xs[j] + q * h) for q in range(grid_per_span + 1)]
        for q in range(grid_per_span):
            u, v = vals[q], vals[q + 1]
            if u >= 0 and v >= 0:
                pos_area += (u + v) * h / 2
            elif u <= 0 and v <= 0:
                neg_area += (u + v) * h / 2
            else:
                r = u / (u - v)
                if u > 0:
                    pos_area += u * r * h / 2
                    neg_area += v * (1 - r) * h / 2
                else:
                    neg_area += u * r * h / 2
                    pos_area += v * (1 - r) * h / 2
        span_min.append(min(vals))
        eta_max = max(eta_max, max(vals))
    negs = sorted(v for v in span_min if v < 0)
    lane_pos = TW_LANE_W * pos_area + TW_LANE_PM * eta_max
    lane_neg = TW_LANE_W * neg_area + TW_LANE_PM * sum(negs[:2])
    pos, neg = max(tp, lane_pos), min(tn, lane_neg)
    return ContLiveMoment(x, tp, tn, lane_pos, lane_neg, pos, neg,
                          taiwan_impact(cont_impact_length(spans, x, +1)),
                          taiwan_impact(cont_impact_length(spans, x, -1)))


def taiwan_cont_live_moment(spans: Sequence[float], x: float, step: float = 0.25,
                            grid_per_span: int = 400) -> ContLiveMoment:
    """斷面 x 的台灣 HS20-44 每車道活載彎矩（不含衝擊）＋對應衝擊係數。"""
    return _live_at(_ILCache(spans), x, step, grid_per_span)


@dataclass
class ContEnvelopeRow:
    x: float
    M_dc: float
    M_dw: float
    M_ll_pos: float     # 含衝擊、車道數與折減
    M_ll_neg: float
    Ms_pos: float       # Service I = DC + DW + LL
    Ms_neg: float
    Mu_pos: float       # Strength I = 1.25DC + 1.50DW + 1.75LL
    Mu_neg: float
    I_pos: float        # 正彎矩衝擊（該跨跨徑）
    I_neg: float        # 負彎矩衝擊（相鄰兩載重跨平均）


def taiwan_cont_envelope(spans: Sequence[float], w_dc: float, w_dw: float, lanes: int,
                         n_per_span: int = 20, step: float = 0.25,
                         grid_per_span: int = 400) -> List[ContEnvelopeRow]:
    """連續梁 DL＋台灣 HS20-44 活載彎矩包絡（不含預力 M2，由呼叫端另加、載重因數 1.0）。"""
    c = _ILCache(spans)
    xs = c.xs
    fac_l = lanes * taiwan_lane_reduction(lanes)
    pts = [0.0]
    for i, L in enumerate(spans):
        for k in range(1, n_per_span + 1):
            pts.append(xs[i] + L * k / n_per_span)
    rows = []
    for x in pts:
        dc, dw = cont_dl_moment(spans, w_dc, x), cont_dl_moment(spans, w_dw, x)
        lv = _live_at(c, x, step, grid_per_span)
        lp, ln = lv.pos * (1 + lv.I_pos) * fac_l, lv.neg * (1 + lv.I_neg) * fac_l
        rows.append(ContEnvelopeRow(x, dc, dw, lp, ln, dc + dw + lp, dc + dw + ln,
                                    1.25 * dc + 1.5 * dw + 1.75 * lp,
                                    1.25 * dc + 1.5 * dw + 1.75 * ln, lv.I_pos, lv.I_neg))
    return rows
