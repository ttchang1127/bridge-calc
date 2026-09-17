"""鋼腱線形設計（G1）：拋物線等效荷載法、荷重平衡率、摩擦損失、曲率半徑。

對應公式卡_鋼腱線形設計、算例_鋼腱線形設計。
平衡荷載法：拋物線鋼腱（垂度 a）對混凝土施加向上等效均佈荷重 w_eq = 8·P·a/L²，
用以抵銷靜載 w_DL；荷重平衡率 LBR = w_eq/w_DL 目標 0.75~0.90（施拉可暫時過平衡）。

單位：力 N、長度 mm（w_eq 回傳 N/mm，數值等於 kN/m）。
摩擦損失介面採規範慣用單位：k 以 /m、路徑長以 m。
"""
import math
from dataclasses import dataclass


def equivalent_load(P: float, a: float, L: float) -> float:
    """拋物線鋼腱等效均佈荷重 w_eq = 8·P·a/L²。P[N] a[mm] L[mm] → N/mm(=kN/m)。"""
    return 8 * P * a / L ** 2


def end_slope(a: float, L: float) -> float:
    """端部傾角 θ_end = 4a/L（對稱拋物線、兩端同高），rad。"""
    return 4 * a / L


def radius_of_curvature(a: float, L: float) -> float:
    """跨中最小曲率半徑 R = L²/(8a)，mm。"""
    return L ** 2 / (8 * a)


def balance_ratio(w_eq: float, w_DL: float) -> float:
    """荷重平衡率 LBR = w_eq / w_DL。"""
    return w_eq / w_DL


def friction_loss(alpha: float, mu: float, k_per_m: float, L_path_m: float) -> float:
    """後張摩擦損失率 ΔP/P = 1 − e^−(μ·α + k·L)。α 累積角變化[rad]、k[/m]、L[m]。"""
    return 1 - math.exp(-(mu * alpha + k_per_m * L_path_m))


@dataclass
class TendonProfileResult:
    a: float               # 垂度 mm
    theta_end: float       # 端部傾角 rad
    R: float               # 曲率半徑 mm
    w_eq_transfer: float   # 施拉時等效荷重 N/mm(=kN/m)
    w_eq_service: float    # 服務時等效荷重 N/mm(=kN/m)
    LBR_transfer: float    # 施拉荷重平衡率
    LBR_service: float     # 服務荷重平衡率
    fric_single_end: float # 單端張拉、遠端（全長）損失率
    fric_dual_mid: float   # 雙端張拉、跨中損失率
    R_ok: bool             # R ≥ R_min


def tendon_profile(Pi: float, Pe: float, a: float, L: float, w_DL: float,
                   mu: float = 0.25, k_per_m: float = 0.003,
                   R_min: float = 3000) -> TendonProfileResult:
    """完整 G1 線形驗算。

    Pi/Pe[N]：施拉/有效預力；a[mm]：跨中垂度（= e_m − e_support，端部偏心 0 時即 e_m）；
    L[mm]：跨徑；w_DL[N/mm]：靜載（簡支可由 8(M_DC+M_DW)/L² 反推）。
    摩擦：單端張拉取全長 α=2θ_end、路徑 L；雙端張拉取跨中 α=θ_end、路徑 L/2。
    """
    th = end_slope(a, L)
    R = radius_of_curvature(a, L)
    w_t = equivalent_load(Pi, a, L)
    w_s = equivalent_load(Pe, a, L)
    L_m = L / 1000.0
    return TendonProfileResult(
        a=a, theta_end=th, R=R,
        w_eq_transfer=w_t, w_eq_service=w_s,
        LBR_transfer=balance_ratio(w_t, w_DL),
        LBR_service=balance_ratio(w_s, w_DL),
        fric_single_end=friction_loss(2 * th, mu, k_per_m, L_m),
        fric_dual_mid=friction_loss(th, mu, k_per_m, L_m / 2.0),
        R_ok=(R >= R_min),
    )


# ── 管道實配排列（G1 §7.1/7.2 構造規定）────────────────────────────
# 分析用的 CGS（鋼腱合力中心）是「把 n 組腱視為一點」的簡化；實際每腹板要排
# n_col 列 × n_row 層的管道，最底層必然低於 CGS。本函式把實配排出來並回推真實
# 的偏心上限 e_max——它一定小於算例以單點算的 ȳb − c − r_duct。
#
# 淨間距需求（台灣 §8.25.2，見公式卡_有效預力與鋼腱配置 §4-1）：
#   s ≥ max(40 mm, 1.5 × 最大骨材粒徑, 孔道外徑)   ← 取大值
# ⚠ 公式卡_鋼腱線形設計 §7.1 表列台灣僅「≥ 40 mm」未含「≥孔道外徑」，兩卡敘述
#   不一致；本引擎取嚴格者（含孔道外徑），偏安全。

@dataclass
class DuctLayoutResult:
    n_per_web: int         # 控制腹板的腱數（= ceil(n/n_web)）
    per_web: list          # 各腹板分配腱數
    n_col: int             # 每層可排列數（腹板內）
    n_row: int             # 層數
    rows: list             # [(y_mm 距梁底, 該層腱數)]，由下往上
    x_offsets: list        # 各列相對腹板中心的橫向偏移 mm
    pitch_v: float         # 層距（= od + s_v）mm
    s_req: float           # 淨間距需求 mm
    s_h: float             # 實際水平淨距 mm（單列時為腹板餘裕）
    s_v: float             # 採用垂直淨距 mm
    web_avail: float       # 腹板可用寬 = web_t − 2·cover
    y_cgs: float           # 鋼腱合力中心距梁底 mm（輸入，實配形心對齊此值）
    y_bot: float           # 最底層腱中心距梁底 mm
    y_top: float           # 最頂層腱中心距梁底 mm
    cover_bot: float       # 最底層管外緣距梁底 mm
    cover_ok: bool
    s_v_ok: bool
    s_h_ok: bool
    top_ok: bool           # 最頂層管外緣仍在斷面內（距頂緣 ≥ cover）
    fits: bool
    e_max: float           # 實配可行的最大偏心 mm（None 表未給 y_b）
    e_max_point: float     # 單點簡化上限 ȳb − cover − od/2 mm（對照用）


def duct_spacing_required(duct_od: float, d_agg: float = 25.0, rule: str = "tw") -> float:
    """管道最小淨間距需求，mm。

    rule="tw"（預設）：max(40, 1.5·d_agg)——台灣橋梁設計規範 §8.25.2 明文的兩個條件，
        亦與公式卡_鋼腱線形設計 §7.1 表列一致。
    rule="od"：再取 max(…, 孔道外徑)——公式卡_有效預力與鋼腱配置 §4-1 的台灣欄另列
        「≥ 孔道外徑」，與 AASHTO「孔道 OD > 100 mm 時淨距 ≥ 孔道外徑」同義。兩卡敘述
        不一致，原規範條文尚未查證，故做成可選；取 "od" 為保守側。

    ⚠ 這個選擇直接決定腹板內排得下幾列：φ100 管、腹板 350、保護層 40 時，
    "tw" → 兩列（淨距需求 40），"od" → 單列（淨距需求 100）。
    """
    base = max(40.0, 1.5 * d_agg)
    return max(base, duct_od) if rule == "od" else base


def duct_layout(n_tendons: int, n_web: int, y_cgs: float,
                duct_od: float = 100.0, web_t: float = 350.0,
                cover: float = 40.0, s_v: float = 40.0,
                d_agg: float = 25.0, y_b: float = None, h: float = None,
                rule: str = "tw") -> DuctLayoutResult:
    """把 n 組鋼腱實際排進腹板，回傳排列幾何與構造檢核。

    n_tendons：鋼腱組數；n_web：腹板數；y_cgs[mm]：鋼腱合力中心距梁底（= ȳb − e）。
    duct_od/web_t/cover/s_v/d_agg[mm]：管道外徑、腹板厚、清保護層、採用垂直淨距、最大骨材粒徑。
    預設 cover=40（台灣 §8.25.1 一般規定；PTBG 實務建議 75）、s_v=40、rule="tw"（見
    duct_spacing_required）。
    y_b[mm]：形心距底（給了才回推 e_max）；h[mm]：梁深（給了才檢核頂側）。

    排列規則：每腹板由下往上逐層填滿，層距 = od + s_v；列數由腹板可用寬與淨間距
    需求決定。排完後整體平移，使實配形心恰等於 y_cgs——引擎用的偏心 e 不因實配改變。
    """
    n_web = max(1, int(n_web))
    n = max(1, int(n_tendons))
    base, rem = divmod(n, n_web)
    per_web = [base + (1 if i < rem else 0) for i in range(n_web)]
    n_per_web = max(per_web)

    s_req = duct_spacing_required(duct_od, d_agg, rule)
    avail = web_t - 2 * cover
    n_col = int((avail + s_req) // (duct_od + s_req)) if avail >= duct_od else 0
    n_col = max(1, n_col)                      # 至少畫一列（放不下時由 s_h_ok 報 ✗）
    s_h = (avail - n_col * duct_od) / (n_col - 1) if n_col > 1 else avail - duct_od
    s_h_ok = (avail >= duct_od) and (n_col == 1 or s_h >= s_req - 1e-9)

    n_row = -(-n_per_web // n_col)             # ceil
    pitch_v = duct_od + s_v
    counts = [min(n_col, n_per_web - i * n_col) for i in range(n_row)]
    rel_cgs = sum(c * i * pitch_v for i, c in enumerate(counts)) / n_per_web
    shift = y_cgs - rel_cgs                    # 平移使實配形心 = 分析用 CGS
    rows = [(i * pitch_v + shift, counts[i]) for i in range(n_row)]

    pitch_h = duct_od + (s_h if n_col > 1 else 0.0)
    x_offsets = [(j - (n_col - 1) / 2.0) * pitch_h for j in range(n_col)]

    y_bot, y_top = rows[0][0], rows[-1][0]
    cover_bot = y_bot - duct_od / 2
    cover_ok = cover_bot >= cover - 1e-9
    s_v_ok = s_v >= s_req - 1e-9
    top_ok = True if h is None else (y_top + duct_od / 2 <= h - cover + 1e-9)
    e_max = None if y_b is None else y_b - (cover + duct_od / 2 + rel_cgs)
    e_max_point = None if y_b is None else y_b - cover - duct_od / 2

    return DuctLayoutResult(
        n_per_web=n_per_web, per_web=per_web, n_col=n_col, n_row=n_row,
        rows=rows, x_offsets=x_offsets, pitch_v=pitch_v,
        s_req=s_req, s_h=s_h, s_v=s_v, web_avail=avail,
        y_cgs=y_cgs, y_bot=y_bot, y_top=y_top, cover_bot=cover_bot,
        cover_ok=cover_ok, s_v_ok=s_v_ok, s_h_ok=s_h_ok, top_ok=top_ok,
        fits=(cover_ok and s_v_ok and s_h_ok and top_ok),
        e_max=e_max, e_max_point=e_max_point)


# ── 摩擦損失沿長度分佈（張拉端配置）────────────────────────────
# ΔP/P = 1 − e^−(μ·α + K·x)：α 是**自張拉端累積到該斷面**的角變化、x 是該段路徑長。
# 兩者都從張拉端起算，所以張拉端配置直接決定每個斷面的損失。
#   α(x) = ∫|e″|dx——反曲時 e″ 變號，必須逐段取絕對值累加，不能用斜率差（那會相消）。
# 配置：'start'/'end' 單端；'both' 雙端（每腱兩端都拉，該點取損失小者）；
#      'alt' 交錯端（半數腱自左、半數自右 → 該斷面取兩者平均）。
# ⚠ 跨中在四種配置下數值相同（左右對稱、路徑各半），差異出現在其餘斷面——
#   單端張拉的遠端損失最大。


@dataclass
class FrictionProfileResult:
    xs: list            # 取樣里程 m
    ratios: list        # 各點摩擦損失率
    at_mid: float       # 跨中
    at_start: float     # 起點端
    at_end: float       # 終點端
    avg: float          # 全長平均（梯形積分）
    max_ratio: float    # 最大損失率
    x_max: float        # 其位置 m
    jack: str


def parabolic_curv_segs(L: float, a: float) -> list:
    """單跨拋物線 e=4a·x(L−x)/L² 的曲率段：κ=|e″|=8a/(L²·1000) rad/m（常數）。

    L[m]、a[mm]（垂度）。回傳 [(x1, x2, kappa)]，供 friction_angle 累加。
    """
    return [(0.0, L, 8 * a / (L * L * 1000.0))]


def friction_angle(segs: list, x: float, L: float, from_end: bool = False) -> float:
    """自張拉端累積到里程 x 的角變化 α（rad）。segs: [(x1, x2, κ[rad/m])]。"""
    tot = 0.0
    for (a, b, k) in segs:
        lo, hi = (max(a, x), min(b, L)) if from_end else (max(a, 0.0), min(b, x))
        if hi > lo:
            tot += k * (hi - lo)
    return tot


def friction_at(x: float, L: float, segs: list, mu: float = 0.25,
                K: float = 0.003, jack: str = "both") -> float:
    """里程 x 的摩擦損失率（0~1）。jack: 'start'|'end'|'both'|'alt'。"""
    lS = 1 - math.exp(-(mu * friction_angle(segs, x, L, False) + K * x))
    lE = 1 - math.exp(-(mu * friction_angle(segs, x, L, True) + K * (L - x)))
    if jack == "start":
        return lS
    if jack == "end":
        return lE
    if jack == "both":
        return min(lS, lE)
    return 0.5 * (lS + lE)                      # alt


def friction_profile(L: float, segs: list, mu: float = 0.25, K: float = 0.003,
                     jack: str = "both", n: int = 20) -> FrictionProfileResult:
    """沿全長取樣的摩擦損失分佈與統計（平均以梯形積分）。"""
    xs = [L * i / n for i in range(n + 1)]
    rs = [friction_at(x, L, segs, mu, K, jack) for x in xs]
    area = sum((rs[i] + rs[i + 1]) / 2 * (xs[i + 1] - xs[i]) for i in range(n))
    imax = max(range(n + 1), key=lambda i: rs[i])
    return FrictionProfileResult(
        xs=xs, ratios=rs,
        at_mid=friction_at(L / 2, L, segs, mu, K, jack),
        at_start=rs[0], at_end=rs[-1], avg=area / L,
        max_ratio=rs[imax], x_max=xs[imax], jack=jack)


# ── 逐腱合成：各腱依自己的張拉端算 Pe，再合成總力與合力位置 ────────
# 為什麼要逐腱：總 Pe 用「平均損失率×總面積」算其實等價（各項對 f_pe 都是線性），
# 真正的差別在**合力位置**——各腱 Pe 不同且高度不同時，合力中心會偏離幾何形心：
#     e_eff = Σ(Pe_i·e_i) / ΣPe_i ≠ e_geom
# 交錯端張拉時，自起點與自終點的腱在同一斷面損失不同；若交錯剛好落在不同層
# （單列多層時序號＝層號），低層與高層的 Pe 就不同 → e_eff 垂直偏移。
# 同層左右配對交錯（多列時）則垂直相消，只剩橫向偏心 x_eff。


@dataclass
class TendonForceResult:
    Pe_total: float     # N
    e_eff: float        # mm，正＝形心下方（以 Pe 加權）
    e_geom: float       # mm，幾何形心（不加權）
    de: float           # e_eff − e_geom（mm）
    x_eff: float        # 橫向合力位置 mm（對稱配置應為 0）
    Pe_avg_ratio: float # 以平均損失率算的總 Pe（對照用，應與 Pe_total 相同）
    per: list           # 每腱 dict：no/y/x/jack/ratio/fpe/Pe


def tendon_forces(tendons: list, x: float, L: float, segs: list, y_b: float,
                  fpj: float, Ap_each: float, other_loss: float,
                  mu: float = 0.25, K: float = 0.003) -> TendonForceResult:
    """逐腱算 Pe 與合力位置。

    tendons：[{"no":str, "y":距梁底 mm, "x":橫向 mm, "jack":'start'|'end'|'both'}]
    x[m]：控制斷面里程；L[m]、segs：線形（見 friction_at）；y_b[mm]：斷面形心距底；
    fpj[MPa]、Ap_each[mm²]（每腱鋼腱面積）、other_loss[MPa]（非摩擦損失合計）。
    """
    per, sP, sPe_e, sPe_x, s_ratio = [], 0.0, 0.0, 0.0, 0.0
    for t in tendons:
        r = friction_at(x, L, segs, mu, K, t.get("jack", "both"))
        fpe = fpj - fpj * r - other_loss
        Pe = fpe * Ap_each
        e_i = y_b - t["y"]
        per.append({"no": t.get("no", ""), "y": t["y"], "x": t.get("x", 0.0),
                    "jack": t.get("jack", "both"), "ratio": r, "fpe": fpe, "Pe": Pe})
        sP += Pe
        sPe_e += Pe * e_i
        sPe_x += Pe * t.get("x", 0.0)
        s_ratio += r
    n = len(tendons)
    e_geom = y_b - sum(t["y"] for t in tendons) / n
    fpe_avg = fpj - fpj * (s_ratio / n) - other_loss
    return TendonForceResult(
        Pe_total=sP, e_eff=(sPe_e / sP if sP else 0.0), e_geom=e_geom,
        de=(sPe_e / sP if sP else 0.0) - e_geom,
        x_eff=(sPe_x / sP if sP else 0.0),
        Pe_avg_ratio=fpe_avg * Ap_each * n, per=per)


def assign_jack(n_per_web: int, mode: str = "alt") -> list:
    """依配置給每腱的張拉端：'alt' 交錯（序號偶數自起點、奇數自終點）／
    'both' 雙端／'start'／'end'。回傳長度 n_per_web 的字串陣列。"""
    if mode in ("both", "start", "end"):
        return [mode] * n_per_web
    return ["start" if i % 2 == 0 else "end" for i in range(n_per_web)]


@dataclass
class TopSlabTendonResult:
    e_hi: float          # 可行偏心上界（形心下為正；管位於頂板最下緣時）
    e_lo: float          # 可行偏心下界（最負；管位於頂板最上緣時）
    y_center: float      # 管中心距梁底 mm
    in_slab: bool        # 管（含保護層）完全位於頂板內
    top_cover: float     # 管外緣距梁頂 mm
    bot_cover: float     # 管外緣距頂板底面 mm
    x_offsets: list      # 各腱水平位置（相對箱梁中心）mm
    s_clear: float       # 水平淨距 mm
    s_req: float
    s_ok: bool
    ok: bool


def top_slab_tendon_check(h: float, yb: float, top_t: float, e_pier: float, n: int,
                          width: float, duct_od: float = 100.0, cover: float = 40.0,
                          d_agg: float = 25.0, rule: str = "tw") -> TopSlabTendonResult:
    """頂板腱（墩頂局部負彎矩腱）於墩頂斷面的構造檢核。

    垂直：管中心 y = ȳb − e 須使管外緣距梁頂 ≥ cover、距頂板底面 ≥ cover
        → e ∈ [ȳb − (h − cover − OD/2), ȳb − (h − top_t + cover + OD/2)]
        例：h 2,100、ȳb 1,329、頂板 250、φ100、保護層 75 → 上下界皆 −646（頂板剛好容納一層）。
    水平：n 束於寬度 width 內等分（間距 width/n、置中），淨距 = 間距 − OD ≥ duct_spacing_required。
        width 由呼叫端決定（建議取兩腹板內緣間，避開腹板內的全長腱）。
    """
    y = yb - e_pier
    e_lo = yb - (h - cover - duct_od / 2)
    e_hi = yb - (h - top_t + cover + duct_od / 2)
    top_cov = h - (y + duct_od / 2)
    bot_cov = (y - duct_od / 2) - (h - top_t)
    in_slab = top_cov >= cover - 1e-6 and bot_cov >= cover - 1e-6
    n = max(1, int(n))
    pitch = width / n
    xo = [(-width / 2 + pitch * (i + 0.5)) for i in range(n)]
    s_clear = pitch - duct_od if n > 1 else float("inf")
    s_req = duct_spacing_required(duct_od, d_agg, rule)
    s_ok = s_clear >= s_req - 1e-6
    return TopSlabTendonResult(e_hi, e_lo, y, in_slab, top_cov, bot_cov, xo,
                               s_clear, s_req, s_ok, in_slab and s_ok)
