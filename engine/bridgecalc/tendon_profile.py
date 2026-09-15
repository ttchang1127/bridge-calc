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
