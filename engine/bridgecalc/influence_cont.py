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


# ── HS20-44 設計卡車：中後軸距 V 可調 4.25～9.15 m（§3.6）──
# 軸序 36 —4.25— 144 —V— 144，雙向行駛。簡支跨恆以 V=4.25 最不利（影響線同號、軸越近越大），
# 故 influence.py 簡支函式維持固定軸距；連續梁兩後軸可能分落兩跨負區，須掃描 V。
TW_REAR_MIN, TW_REAR_MAX, TW_REAR_STEP = 4.25, 9.15, 0.25
_GRID = 0.05                                     # 影響線取樣格距（m）；軸距 4.25、9.15 皆為其整數倍


def taiwan_rear_spacings(step: float = TW_REAR_STEP) -> List[float]:
    out, v = [], TW_REAR_MIN
    while v < TW_REAR_MAX - 1e-9:
        out.append(round(v, 6))
        v += step
    out.append(TW_REAR_MAX)
    return out


def _truck_sets(V: float):
    P = TW_HS20_AXLES
    d = (0.0, 4.25, 4.25 + V)
    return ((P, d), (tuple(reversed(P)), tuple(d[-1] - t for t in reversed(d))))


def _truck_extremes(vals: Sequence[float], tot: float, at_x=None, s_every: int = 2,
                    spacings: Sequence[float] = None):
    """vals[k] = η(k·_GRID)；回傳 (最大, 最小, 最大時 V, 最小時 V)。

    at_x：(x, eta_side_fn) 時另將各軸逐一放在斷面 x 上（左右極限）評估，
    補足格點未必落在斷面上（剪力影響線於斷面跳 1、彎矩峰值在斷面）的情形。
    """
    N = len(vals) - 1
    best_p = best_n = 0.0
    vp = vn = TW_REAR_MIN
    for V in (spacings or taiwan_rear_spacings()):
        for P_set, d_set in _truck_sets(V):
            di = [int(round(t / _GRID)) for t in d_set]
            span_i = di[-1]
            for k in range(-span_i, N + 1, s_every):
                v = 0.0
                for P, o in zip(P_set, di):
                    j = k + o
                    if 0 <= j <= N:
                        v += P * vals[j]
                if v > best_p:
                    best_p, vp = v, V
                if v < best_n:
                    best_n, vn = v, V
            if at_x is not None:
                x, fn = at_x
                for a_idx in range(3):
                    s0 = x - d_set[a_idx]
                    for side in (-1, 1):
                        v = 0.0
                        for i2, (P, t) in enumerate(zip(P_set, d_set)):
                            ax = s0 + t
                            if -1e-9 <= ax <= tot + 1e-9:
                                v += P * fn(min(max(ax, 0.0), tot), side if i2 == a_idx else 0)
                        if v > best_p:
                            best_p, vp = v, V
                        if v < best_n:
                            best_n, vn = v, V
    return best_p, best_n, vp, vn


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
    V_pos: float = 4.25  # 卡車控制時之中後軸距
    V_neg: float = 4.25


class _ILCache:
    def __init__(self, spans, stiff=None):
        self.spans = list(spans)
        self.xs = _supports(spans)
        self.tot = self.xs[-1]
        self._X = {}
        self.flex = None
        if stiff is not None:
            from .variable_section import ContFlex
            self.flex = ContFlex(spans, stiff.I_rel, stiff.breaks)

    def X(self, p):
        key = round(p, 6)
        v = self._X.get(key)
        if v is None:
            v = self._X[key] = (self.flex.support_moments_point(p) if self.flex
                                else cont_support_moments_point(self.spans, p))
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
    # 設計卡車：雙向、中後軸距 V 掃描；影響線取 0.05 m 格點，另將各軸放在斷面上補峰值
    N = int(round(tot / _GRID))
    vals = [c.eta(x, min(k * _GRID, tot)) for k in range(N + 1)]
    tp, tn, Vp_, Vn_ = _truck_extremes(vals, tot, (x, lambda p, sd: c.eta(x, p)))
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
                          taiwan_impact(cont_impact_length(spans, x, -1)), Vp_, Vn_)


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
                         grid_per_span: int = 400, stiff=None) -> List[ContEnvelopeRow]:
    """連續梁 DL＋台灣 HS20-44 活載彎矩包絡（不含預力 M2，由呼叫端另加、載重因數 1.0）。

    stiff：變斷面剖面（variable_section.haunch_profile）→ 影響線以變 EI 柔度求解、
    自重 w_dc 依斷面積比 A(x)/A_ref 放大（附加恆載 w_dw 不變）。未給時為等斷面閉合解。
    """
    c = _ILCache(spans, stiff)
    xs = c.xs
    fac_l = lanes * taiwan_lane_reduction(lanes)
    pts = [0.0]
    for i, L in enumerate(spans):
        for k in range(1, n_per_span + 1):
            pts.append(xs[i] + L * k / n_per_span)
    if c.flex is not None:
        wdc_fn = lambda t: w_dc * stiff.A_rel(t)
        wdw_fn = lambda t: w_dw
        Xdc, Xdw = c.flex.support_moments_dist(wdc_fn), c.flex.support_moments_dist(wdw_fn)
    rows = []
    for x in pts:
        if c.flex is None:
            dc, dw = cont_dl_moment(spans, w_dc, x), cont_dl_moment(spans, w_dw, x)
        else:
            dc = c.flex.moment_at(Xdc, x, c.flex.M0_dist(wdc_fn, x))
            dw = c.flex.moment_at(Xdw, x, c.flex.M0_dist(wdw_fn, x))
        lv = _live_at(c, x, step, grid_per_span)
        lp, ln = lv.pos * (1 + lv.I_pos) * fac_l, lv.neg * (1 + lv.I_neg) * fac_l
        rows.append(ContEnvelopeRow(x, dc, dw, lp, ln, dc + dw + lp, dc + dw + ln,
                                    1.25 * dc + 1.5 * dw + 1.75 * lp,
                                    1.25 * dc + 1.5 * dw + 1.75 * ln, lv.I_pos, lv.I_neg))
    return rows


# ════════════════ 連續梁剪力影響線與包絡 ════════════════
# 斷面 x 的剪力 V = V0（簡支，載重與斷面同跨）＋(X_{i+1} − X_i)/L_i（支承彎矩梯度）。
# 正負號同 dM/dx（左段向上為正）。斷面處剪力影響線跳 1：side="L"／"R" 指斷面取左側或右側極限；
# 恰在內支承時 "L" 屬左跨末端、"R" 屬右跨起點（兩者相差該支承反力）。
# 台灣：車道集中載重（剪力）116 kN、一個；衝擊長度＝該點至所屬跨較遠支點之距離（§3.13）。

def _shear_span(xs, spans, x, side):
    n = len(spans)
    for j in range(1, n):
        if abs(x - xs[j]) < 1e-9:
            return j - 1 if side == "L" else j
    return _span_of(xs, x)


def _V0_point(xs, spans, i, x, p, p_side=0):
    """簡支跨 i 於斷面 x 的剪力影響值；p==x 時 p_side −1＝載重在斷面左、+1＝在右。"""
    L = spans[i]
    a, q = x - xs[i], p - xs[i]
    if q < -1e-9 or q > L + 1e-9:
        return 0.0
    if abs(q - a) < 1e-9:
        return -q / L if p_side < 0 else (L - q) / L
    return -q / L if q < a else (L - q) / L


def cont_shear_il(spans: Sequence[float], x: float, p: float, side: str = "R",
                  p_side: int = 0, _X=None) -> float:
    """斷面 x（side 側）的剪力影響線縱距（單位載重於 p）。"""
    xs = _supports(spans)
    i = _shear_span(xs, spans, x, side)
    X = _X if _X is not None else cont_support_moments_point(spans, p)
    in_span = xs[i] - 1e-9 <= p <= xs[i + 1] + 1e-9
    V0 = _V0_point(xs, spans, i, x, p, p_side if p_side else (1 if side == "R" else -1)) if in_span else 0.0
    return V0 + (X[i + 1] - X[i]) / spans[i]


def cont_dl_shear(spans: Sequence[float], w: float, x: float, side: str = "R", flex=None,
                  w_fn=None, X=None) -> float:
    """均布（或分佈 w_fn）恆載於斷面 x（side 側）的剪力。"""
    xs = _supports(spans)
    i = _shear_span(xs, spans, x, side)
    L, a = spans[i], x - xs[i]
    if flex is None:
        Xs = cont_support_moments_uniform(spans, [w] * len(spans))
        V0 = w * (L / 2 - a)
    else:
        from .variable_section import gauss_integrate
        wf = w_fn or (lambda t: w)
        Xs = X if X is not None else flex.support_moments_dist(wf)
        RA = gauss_integrate(lambda t: wf(t) * (L - (t - xs[i])) / L, xs[i], xs[i + 1], flex.breaks, 0.25)
        V0 = RA - gauss_integrate(wf, xs[i], x, flex.breaks, 0.25)
    return V0 + (Xs[i + 1] - Xs[i]) / L


def cont_shear_impact_length(spans: Sequence[float], x: float, side: str = "R") -> float:
    """剪力衝擊長度：該點至所屬跨較遠支點之距離。"""
    xs = _supports(spans)
    i = _shear_span(xs, spans, x, side)
    a = x - xs[i]
    return max(a, spans[i] - a)


@dataclass
class ContLiveShear:
    x: float
    side: str
    truck_pos: float
    truck_neg: float
    lane_pos: float
    lane_neg: float
    pos: float
    neg: float
    I: float


def _live_shear_at(c: "_ILCache", x: float, side: str, step: float, grid_per_span: int) -> ContLiveShear:
    from .influence import TW_LANE_PV
    spans, xs, tot = c.spans, c.xs, c.tot
    i = _shear_span(xs, spans, x, side)

    def eta(p, ps=0):
        return cont_shear_il(spans, x, p, side, ps, c.X(p))
    N = int(round(tot / _GRID))
    vals = [eta(min(k * _GRID, tot), (1 if k * _GRID > x else -1) if abs(k * _GRID - x) < 1e-9 else 0)
            for k in range(N + 1)]
    tp, tn, _, _ = _truck_extremes(vals, tot, (x, lambda p, sd: eta(p, sd)))
    pos_area = neg_area = 0.0
    eta_max = eta_min = 0.0
    for j, L in enumerate(spans):
        a0, b0 = xs[j], xs[j + 1]
        pieces = [(a0, b0)]
        if j == i and a0 < x < b0:
            pieces = [(a0, x), (x, b0)]
        for u0, v0 in pieces:
            n = max(1, int(round(grid_per_span * (v0 - u0) / L)))
            h = (v0 - u0) / n
            vals = []
            for q in range(n + 1):
                p = u0 + q * h
                ps = -1 if (j == i and abs(p - x) < 1e-9 and q == n) else (1 if (j == i and abs(p - x) < 1e-9) else 0)
                vals.append(eta(p, ps))
            for q in range(n):
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
            eta_max, eta_min = max(eta_max, max(vals)), min(eta_min, min(vals))
    lane_pos = TW_LANE_W * pos_area + TW_LANE_PV * eta_max
    lane_neg = TW_LANE_W * neg_area + TW_LANE_PV * eta_min
    return ContLiveShear(x, side, tp, tn, lane_pos, lane_neg, max(tp, lane_pos), min(tn, lane_neg),
                         taiwan_impact(cont_shear_impact_length(spans, x, side)))


def taiwan_cont_live_shear(spans: Sequence[float], x: float, side: str = "R", step: float = 0.25,
                           grid_per_span: int = 400, stiff=None) -> ContLiveShear:
    """斷面 x（side 側）的台灣 HS20-44 每車道活載剪力（不含衝擊）＋衝擊係數。"""
    return _live_shear_at(_ILCache(spans, stiff), x, side, step, grid_per_span)


@dataclass
class ContShearRow:
    x: float
    side: str
    V_dc: float
    V_dw: float
    V_ll_pos: float     # 含衝擊、車道數與折減
    V_ll_neg: float
    Vu_pos: float       # 1.25DC + 1.50DW + 1.75LL（不含預力次剪力 V2）
    Vu_neg: float
    I: float


def taiwan_cont_shear_at(spans: Sequence[float], x: float, side: str, w_dc: float, w_dw: float,
                         lanes: int, step: float = 0.25, grid_per_span: int = 400, stiff=None,
                         _cache=None) -> ContShearRow:
    """連續梁斷面 x（side 側）DL＋台灣 HS20-44 活載設計剪力。stiff 同 taiwan_cont_envelope。"""
    c = _cache or _ILCache(spans, stiff)
    fac = lanes * taiwan_lane_reduction(lanes)
    if c.flex is None:
        dc = cont_dl_shear(spans, w_dc, x, side)
        dw = cont_dl_shear(spans, w_dw, x, side)
    else:
        wdc_fn = lambda t: w_dc * stiff.A_rel(t)
        dc = cont_dl_shear(spans, w_dc, x, side, c.flex, wdc_fn)
        dw = cont_dl_shear(spans, w_dw, x, side, c.flex)
    lv = _live_shear_at(c, x, side, step, grid_per_span)
    lp, ln = lv.pos * (1 + lv.I) * fac, lv.neg * (1 + lv.I) * fac
    return ContShearRow(x, side, dc, dw, lp, ln, 1.25 * dc + 1.5 * dw + 1.75 * lp,
                        1.25 * dc + 1.5 * dw + 1.75 * ln, lv.I)


def secondary_shear(fm, spans: Sequence[float], x: float, side: str = "R") -> float:
    """預力次剪力 V2 = dM2/dx = (X_{i+1} − X_i)/L_i（fm：ForceMethodM2Result）。"""
    xs = _supports(spans)
    i = _shear_span(xs, spans, x, side)
    return (fm.X[i + 1] - fm.X[i]) / spans[i]


def design_shear_with_V2(Vu_pos: float, Vu_neg: float, V2: float) -> float:
    """設計剪力（帶號）：Vu = max(|V_載重|, |V_載重 + V2|)（次剪力有利時不折減；算例_端跨 propped cantilever）。"""
    cand = [Vu_pos, Vu_neg, Vu_pos + V2, Vu_neg + V2]
    return max(cand, key=abs)
