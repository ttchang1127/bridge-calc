"""鋼腱線形設計（G1）：拋物線等效荷載法、荷重平衡率、摩擦損失、曲率半徑。

對應公式卡_鋼腱線形設計、算例_鋼腱線形設計。
平衡荷載法：拋物線鋼腱（垂度 a）對混凝土施加向上等效均佈荷重 w_eq = 8·P·a/L²，
用以抵銷靜載 w_DL；荷重平衡率 LBR = w_eq/w_DL 目標 0.75~0.90（施拉可暫時過平衡）。

單位：力 N、長度 mm（w_eq 回傳 N/mm，數值等於 kN/m）。
摩擦損失介面採規範慣用單位：k 以 /m、路徑長以 m。
"""
import math
from typing import Sequence
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
# 淨間距需求（台灣 §8.25.2，2026-09-18 已核原文）：
#   s ≥ max(40 mm, 1.5 × 最大骨材粒徑)   ← 明文僅此兩條件，無「≥ 孔道外徑」
# 「≥ 孔道外徑」是 AASHTO 條文（OD > 100 mm），曾誤植入 B3 §4-1 台灣欄，已更正；
# 引擎以 rule="od" 保留為非規範要求的保守設計慣例（見 duct_spacing_required）。
# 另 §8.25.3：彎折／撓曲段套管最多三束可捆紮，但距構材端部 90 cm 內須維持上列間距。

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
    cover_side: float = None  # 最外側管壁至腹板側面的淨距 mm（side_margin=0 時＝cover）


def duct_spacing_required(duct_od: float, d_agg: float = 25.0, rule: str = "tw") -> float:
    """管道最小淨間距需求，mm。

    rule="tw"（預設）：max(40, 1.5·d_agg)——台灣橋梁設計規範 §8.25.2 明文僅此兩條件
        （2026-09-18 核對第八章原文 p.178–179，全章無「≥ 孔道外徑」之規定）。
    rule="od"：再取 max(…, 孔道外徑)——AASHTO「孔道 OD > 100 mm 時淨距 ≥ 孔道外徑」，
        **非台灣規範要求**，保留為欲與 AASHTO 一致或採嚴格品管時的保守設計慣例。

    ⚠ 這個選擇直接決定腹板內排得下幾列：φ100 管、腹板 350、保護層 40 時，
    "tw" → 兩列（淨距需求 40），"od" → 單列（淨距需求 100）。
    """
    base = max(40.0, 1.5 * d_agg)
    return max(base, duct_od) if rule == "od" else base


def duct_layout(n_tendons: int, n_web: int, y_cgs: float,
                duct_od: float = 100.0, web_t: float = 350.0,
                cover: float = 40.0, s_v: float = 40.0,
                d_agg: float = 25.0, y_b: float = None, h: float = None,
                rule: str = "tw", side_margin: float = 0.0) -> DuctLayoutResult:
    """把 n 組鋼腱實際排進腹板，回傳排列幾何與構造檢核。

    n_tendons：鋼腱組數；n_web：腹板數；y_cgs[mm]：鋼腱合力中心距梁底（= ȳb − e）。
    duct_od/web_t/cover/s_v/d_agg[mm]：管道外徑、腹板厚、清保護層、採用垂直淨距、最大骨材粒徑。
    預設 cover=40（台灣 §8.25.1 一般規定；PTBG 實務建議 75）、s_v=40、rule="tw"（見
    duct_spacing_required）。
    y_b[mm]：形心距底（給了才回推 e_max）；h[mm]：梁深（給了才檢核頂側）。

    排列規則：每腹板由下往上逐層填滿，層距 = od + s_v；列數由腹板可用寬與淨間距
    需求決定。排完後整體平移，使實配形心恰等於 y_cgs——引擎用的偏心 e 不因實配改變。

    ⚠ 多列時剩餘寬度全數分給列間距 → 最外側管的側向保護層**恰等於 cover、餘裕為零**。
    side_margin[mm]：兩側各再預留的餘裕（施工誤差、箍筋位置），可用寬 = web_t − 2(cover+side_margin)。
    預設 0＝原行為。實際側向淨保護層見結果 cover_side。
    """
    n_web = max(1, int(n_web))
    n = max(1, int(n_tendons))
    base, rem = divmod(n, n_web)
    per_web = [base + (1 if i < rem else 0) for i in range(n_web)]
    n_per_web = max(per_web)

    s_req = duct_spacing_required(duct_od, d_agg, rule)
    avail = web_t - 2 * (cover + side_margin)
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
        e_max=e_max, e_max_point=e_max_point,
        cover_side=(web_t - (x_offsets[-1] - x_offsets[0]) - duct_od) / 2.0)


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


def tendon_slip_loss(x: float, L: float, segs: list, fpj: float, jack: str = "both",
                     mu: float = 0.25, K: float = 0.003, slip_mm: float = 6.0,
                     Ep: float = 195000.0) -> float:
    """單腱於里程 x 的錨具滑移損失 MPa（線形化摩擦梯度，見 anchor_slip_loss）。

    jack='both' → 兩端各影響半長（L_m=L/2、d＝距最近端）；'start'/'end' → L_m=L、d 自該端量。
    p＝**該滑移端**至 L_m 的摩擦損失 ÷ L_m——雙端張拉時依 x 靠近哪一端取該端的梯度；
    線形不對稱時兩端梯度不同，取固定一端會使段內左右不對稱。slip_mm=0 即回傳 0（預設行為不變）。
    """
    if slip_mm <= 0 or L <= 0:
        return 0.0
    if jack == "both":
        L_m = L / 2
        if x <= L_m:
            d, r = x, friction_at(L_m, L, segs, mu, K, "start")
        else:
            d, r = L - x, friction_at(L_m, L, segs, mu, K, "end")
    elif jack == "end":
        L_m, d = L, L - x
        r = friction_at(0.0, L, segs, mu, K, "end")
    else:
        L_m, d = L, x
        r = friction_at(L, L, segs, mu, K, "start")
    p = fpj * r / L_m if L_m > 0 else 0.0
    return anchor_slip_loss(d, L_m, p, slip_mm, Ep).dsigma


def tendon_forces(tendons: list, x: float, L: float, segs: list, y_b: float,
                  fpj: float, Ap_each: float, other_loss: float,
                  mu: float = 0.25, K: float = 0.003, slip_mm: float = 0.0,
                  Ep: float = 195000.0, anchors: Sequence[float] = None) -> TendonForceResult:
    """逐腱算 Pe 與合力位置。

    tendons：[{"no":str, "y":距梁底 mm, "x":橫向 mm, "jack":'start'|'end'|'both'}]
    x[m]：控制斷面里程；L[m]、segs：線形（見 friction_at）；y_b[mm]：斷面形心距底；
    fpj[MPa]、Ap_each[mm²]（每腱鋼腱面積）、other_loss[MPa]（非摩擦損失合計）。
    slip_mm：錨具滑移量（預設 0＝不計，與既有結果相同）；各腱依自身張拉端計算（見 tendon_slip_loss）。
    anchors：中間錨碇段的錨碇里程（含兩端，如 [0, 40, 80]）——摩擦與滑移只在各段內累積
    （見 segmented_tendon_force）；不給則整條腱為一段。
    """
    per, sP, sPe_e, sPe_x, s_ratio = [], 0.0, 0.0, 0.0, 0.0
    for t in tendons:
        jk = t.get("jack", "both")
        if anchors:
            f = segmented_tendon_force(x, anchors, segs, fpj, Ap_each, other_loss, mu, K, jk, slip_mm, Ep)
            r, sl, fpe = f.friction / fpj, f.slip, f.fpe
        else:
            r = friction_at(x, L, segs, mu, K, jk)
            sl = tendon_slip_loss(x, L, segs, fpj, jk, mu, K, slip_mm, Ep)
            fpe = fpj - fpj * r - sl - other_loss
        Pe = fpe * Ap_each
        e_i = y_b - t["y"]
        per.append({"no": t.get("no", ""), "y": t["y"], "x": t.get("x", 0.0),
                    "jack": jk, "ratio": r, "slip": sl, "fpe": fpe, "Pe": Pe})
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


@dataclass
class AnchorSlipResult:
    dsigma: float        # 該點滑移損失 MPa
    L_set: float         # 影響長度 m（受 L_m 限制時為 L_m）
    dsigma_anchor: float # 錨具端最大損失 MPa
    p: float             # 摩擦損失梯度 MPa/m
    capped: bool         # L_set 超過 L_m（損失圖整體平移）


def anchor_slip_loss(d: float, L_m: float, p: float, slip_mm: float = 6.0,
                     Ep: float = 195000.0) -> AnchorSlipResult:
    """錨具滑移損失（公式卡_後張預力短期損失 §二，簡化法：摩擦損失沿長度線性）。

    d：距張拉端 m；L_m：滑移可影響的最遠距離 m（單端張拉＝全長、雙端張拉＝半長，即摩擦損失最低點）；
    p：摩擦損失梯度 MPa/m。面積條件 ∫Δσ dx = Δa·Ep（mm·MPa → 除 1000 化為 MPa·m）：
      L_set = √(Δa·Ep/(1000·p))；Δσ(d) = 2p·(L_set − d)（d < L_set）。
    L_set > L_m 時損失圖整體平移：Δσ(d) = 2p·(L_m − d) + c，c = (Δa·Ep/1000 − p·L_m²)/L_m。
    """
    A = slip_mm * Ep / 1000.0
    if p <= 0:
        c = A / L_m if L_m > 0 else 0.0
        return AnchorSlipResult(c, L_m, c, p, True)
    Ls = math.sqrt(A / p)
    if Ls <= L_m:
        ds = 2 * p * (Ls - d) if d < Ls else 0.0
        return AnchorSlipResult(ds, Ls, 2 * p * Ls, p, False)
    c = (A - p * L_m ** 2) / L_m
    return AnchorSlipResult(2 * p * (L_m - min(d, L_m)) + c, L_m, 2 * p * L_m + c, p, True)


@dataclass
class CapTendonForce:
    P: float             # 該點總預力 kN（不在腱範圍內為 0）
    fpe: float           # MPa
    friction: float      # 摩擦損失 MPa
    slip: float          # 錨具滑移損失 MPa
    other: float         # 其他（ES、潛變、乾縮、鬆弛）MPa
    s: float             # 距左錨 m
    L_t: float           # 該組腱全長 m


def pier_cap_tendon_force(pc, x: float, fpj: float, Aps: float, other_loss: float,
                          mu: float = 0.25, K: float = 0.003, slip_mm: float = 6.0,
                          Ep: float = 195000.0) -> CapTendonForce:
    """墩頂局部腱（雙端張拉）逐點有效預力：摩擦（依自身線形）＋錨具滑移＋其他損失。

    pc：pier_cap_tendon_segs 結果（每墩兩段）；Aps：該組全部鋼腱面積 mm²。
    摩擦由 friction_at（雙端張拉取兩端損失小者）；滑移依 anchor_slip_loss，L_m＝半長、
    p＝半長處摩擦損失 ÷ 半長（線性化）。other_loss 通常沿用全長腱的非摩擦損失（f_cgp 差異略去）。
    """
    segs = pc.segs
    for k in range(0, len(segs), 2):
        sl, sr = segs[k], segs[k + 1]
        xa, xb = sl.x1, sr.x2
        if xa - 1e-9 <= x <= xb + 1e-9:
            L_t = xb - xa
            local = [(0.0, sl.x2 - xa, 2 * abs(sl.c) / 1000.0), (sl.x2 - xa, L_t, 2 * abs(sr.c) / 1000.0)]
            s = min(max(x - xa, 0.0), L_t)
            r = friction_at(s, L_t, local, mu, K, "both")
            L_m = L_t / 2
            r_mid = friction_at(L_m, L_t, local, mu, K, "start" if s <= L_m else "end")
            p = fpj * r_mid / L_m if L_m > 0 else 0.0
            sl_ = anchor_slip_loss(min(s, L_t - s), L_m, p, slip_mm, Ep).dsigma
            fr = fpj * r
            fpe = fpj - fr - sl_ - other_loss
            return CapTendonForce(fpe * Aps / 1000.0, fpe, fr, sl_, other_loss, s, L_t)
    return CapTendonForce(0.0, 0.0, 0.0, 0.0, other_loss, 0.0, 0.0)


@dataclass
class SegmentedTendonForce:
    P: float             # kN
    fpe: float           # MPa
    friction: float      # MPa
    slip: float          # MPa
    seg: int             # 所在段（0 起）
    L_seg: float         # 該段長度 m
    s: float             # 距該段左端 m


def segmented_tendon_force(x: float, anchors: Sequence[float], curv_segs: Sequence,
                           fpj: float, Aps: float, other_loss: float,
                           mu: float = 0.25, K: float = 0.003, jack: str = "both",
                           slip_mm: float = 0.0, Ep: float = 195000.0) -> SegmentedTendonForce:
    """分段（中間錨碇）鋼腱於里程 x 的有效預力。

    anchors：錨碇里程（含兩端），例如 40+40 於墩頂設中間錨碇 → [0, 40, 80]；
    curv_segs：全域曲率段 [(x1, x2, |e″| rad/m), ...]（同 friction_at）。
    每段各自張拉（jack 對該段而言：'both'／'start'／'end'），摩擦與滑移都只在該段內累積——
    **這就是中間錨碇段降低摩擦損失的原因**：α 與 x 都自該段的張拉端起算（G1 不合格決策樹：>60 m 建議）。
    """
    a = list(anchors)
    for i in range(len(a) - 1):
        if a[i] - 1e-9 <= x <= a[i + 1] + 1e-9:
            x0, x1 = a[i], a[i + 1]
            L = x1 - x0
            loc = []
            for (c0, c1, k) in curv_segs:
                lo, hi = max(c0, x0), min(c1, x1)
                if hi > lo:
                    loc.append((lo - x0, hi - x0, k))
            s = min(max(x - x0, 0.0), L)
            r = friction_at(s, L, loc, mu, K, jack)
            sl = tendon_slip_loss(s, L, loc, fpj, jack, mu, K, slip_mm, Ep)
            fpe = fpj - fpj * r - sl - other_loss
            return SegmentedTendonForce(fpe * Aps / 1000.0, fpe, fpj * r, sl, i, L, s)
    return SegmentedTendonForce(0.0, 0.0, 0.0, 0.0, -1, 0.0, 0.0)


def segmented_friction_profile(anchors: Sequence[float], curv_segs: Sequence, fpj: float,
                               mu: float = 0.25, K: float = 0.003, jack: str = "both",
                               slip_mm: float = 0.0, n: int = 40) -> dict:
    """分段腱沿全長的損失率分佈：回傳 {xs, ratios（摩擦＋滑移佔 fpj）, max_ratio, x_max, avg}。"""
    tot = anchors[-1] - anchors[0]
    xs, rs = [], []
    for i in range(n + 1):
        x = anchors[0] + tot * i / n
        f = segmented_tendon_force(x, anchors, curv_segs, fpj, 1.0, 0.0, mu, K, jack, slip_mm)
        xs.append(x)
        rs.append((f.friction + f.slip) / fpj)
    mx = max(rs)
    return {"xs": xs, "ratios": rs, "max_ratio": mx, "x_max": xs[rs.index(mx)],
            "avg": sum(rs) / len(rs)}


# ── 套管尺寸檢核（台灣 §8.25.4、表 8.3；PTI Manual 6th Table 4.3/4.4）──────
# 🔴 表 8.3 與 PTI Table 4.4 列的都是**內徑**；外徑／壁厚依廠商波紋形狀與材質而定
#   （PTI §4.4.5），規範不提供。排列、淨間距、保護層用**外徑**；面積比用**內徑**——
#   兩者混用會使排列檢核偏不保守（把內徑當外徑）。
TW_DUCT_MAX_ID = {   # 表 8.3 套管最大內徑參考值 mm：{鋼絞線徑: {股數: 內徑}}
    12.7: {22: 90.0, 19: 90.0, 12: 75.0, 7: 55.0},
    15.2: {22: 110.0, 19: 100.0, 12: 85.0, 7: 70.0},
}
DUCT_AREA_RATIO = {  # 套管內面積 / 鋼腱淨面積 下限
    "tw": 2.0,           # 台灣 §8.25.4(1)：多股鋼腱「至少應為鋼腱淨面積之二倍」
    "pti_push": 2.25,    # PTI Table 4.3：strand-push-through
    "pti_pull": 2.5,     # PTI Table 4.3：strand-pull-through
    "pti_short": 2.0,    # PTI：≤ 30 m 短腱或空間受限
}


@dataclass
class DuctSizeResult:
    A_ps: float          # 鋼腱淨面積 mm²
    A_duct: float        # 套管內面積 mm²
    ratio: float         # A_duct / A_ps
    ratio_req: float
    area_ok: bool
    id_max_tw: float     # 表 8.3 最大內徑（查無則 None）
    id_ok: bool          # 未超過表 8.3（查無則 True）
    od_gt_id: bool       # 外徑 > 內徑（否則代表把內徑當外徑用）
    wall: float          # (OD − ID)/2 mm


def duct_size_check(n_strands: int, duct_id: float, duct_od: float = None,
                    strand_dia: float = 15.2, strand_area: float = 140.0,
                    rule: str = "tw") -> DuctSizeResult:
    """套管尺寸檢核：面積比（§8.25.4／PTI 4.3）、表 8.3 最大內徑、外徑是否大於內徑。

    duct_id：內徑 mm（面積比用）；duct_od：外徑 mm（排列用；不給則視為未知）。
    """
    Aps = n_strands * strand_area
    Ad = math.pi * duct_id ** 2 / 4.0
    req = DUCT_AREA_RATIO[rule]
    idm = TW_DUCT_MAX_ID.get(strand_dia, {}).get(n_strands)
    od = duct_id if duct_od is None else duct_od
    return DuctSizeResult(
        A_ps=Aps, A_duct=Ad, ratio=Ad / Aps, ratio_req=req, area_ok=Ad / Aps >= req - 1e-9,
        id_max_tw=idm, id_ok=(idm is None or duct_id <= idm + 1e-9),
        od_gt_id=od > duct_id + 1e-9, wall=(od - duct_id) / 2.0)
