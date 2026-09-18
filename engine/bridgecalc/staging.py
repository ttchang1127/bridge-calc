"""施工階段體系轉換的潛變重分配（Creep Redistribution on Change of Statical System）。

★ 這是「逐跨施工」「先簡支後連續」與「懸臂工法合龍」共同的核心問題，
   而且**實務上常是控制因素**——多數手算把連續梁的恆載直接用連續體系算，
   那只有在「整座橋一次澆築成形」時才對。

物理：載重在時刻 t₀ 作用於**施工當時的**靜力體系 I（逐跨施工＝各跨簡支）；
     到 t₁ 建立連續性，體系變成 II。此後混凝土潛變想繼續依體系 I 的曲率變形，
     但新接頭不准相對轉角發生 → 產生**束制彎矩**，內力朝「載重一開始就作用在
     體系 II」的解漂移：

     M(∞) = M_I + λ·(M_II − M_I)

     λ＝重分配係數，由 t₁ 之後**尚未發生的**潛變量 Δφ = φ(∞,t₀) − φ(t₁,t₀) 決定。

     ⚠ 決定 λ 的是 Δφ 而不是 φ(∞)——**連續性建立得越早，重分配越多**。
       這與直覺相反的地方在於：早連續 ⇒ 剩餘潛變多 ⇒ 墩頂負彎矩大。

     ★ ΔM = λ(M_II − M_I) 是**束制彎矩**，不對應任何外載 → 沿梁**在支承間為線性**。
       這是本模組最強的自我驗證（見 `redistribution_is_linear`）。

單位：彎矩 kN·m、長度 m。
"""
import math
from dataclasses import dataclass
from typing import Sequence


# ── 重分配係數 λ ───────────────────────────────────────────────
@dataclass
class RedistFactor:
    lam: float          # 重分配係數
    dphi: float         # Δφ＝連續後尚未發生的潛變
    chi: float          # 老化係數（Trost 用）
    method: str
    valid: bool         # λ 是否落在 [0,1]（超出即模型適用範圍外）
    note: str


def redistribution_factor(dphi: float, chi: float = 0.8,
                          method: str = "trost") -> RedistFactor:
    """體系轉換的潛變重分配係數。

    method="trost"（預設，AAEM 老化係數法，較準）：

        λ = Δφ / (1 + χ·Δφ)

    由相容條件導得：連續後的潛變轉角 θ₁·Δφ 必須被束制彎矩抵銷，而束制彎矩是
    **逐漸**施加的，故其效應以老化調整模數 E_c/(1+χΔφ) 計算。

    method="dischinger"（rate-of-creep，經典但忽略延遲彈性）：

        λ = 1 − e^(−Δφ)

    兩式為**同一物理的不同近似**，小 Δφ 時一階一致（Δφ=0.1 相對差 2.6%）；
    實務範圍 Δφ ≤ 3 相差約 0.07（Dischinger 偏大）。

    ⚠ Trost 式在 Δφ > 1/(1−χ)（χ=0.8 時為 5.0）會給出 λ>1，即「超過最終體系的值」，
      那已超出模型適用範圍——此時 `valid=False`，應改用數值逐步積分法。
    """
    if dphi < 0:
        raise ValueError("Δφ 不可為負（連續後尚未發生的潛變量）")
    if method == "dischinger":
        lam = 1 - math.exp(-dphi)
    elif method == "trost":
        lam = dphi / (1 + chi * dphi) if (1 + chi * dphi) != 0 else 0.0
    else:
        raise ValueError("method 須為 'trost' 或 'dischinger'")
    valid = 0.0 <= lam <= 1.0 + 1e-12
    note = "" if valid else f"λ={lam:.3f} 超出 [0,1]，Δφ 過大已逸出模型適用範圍"
    return RedistFactor(lam=lam, dphi=dphi, chi=chi, method=method,
                        valid=valid, note=note)


# ── 重分配計算 ────────────────────────────────────────────────
@dataclass
class RedistPoint:
    x: float        # 里程 m
    M_I: float      # 施工體系（as-built）彎矩 kN·m
    M_II: float     # 最終體系（as if built monolithic）彎矩 kN·m
    dM: float       # 束制彎矩 λ(M_II − M_I)
    M_inf: float    # 重分配後 M_I + dM


@dataclass
class StagingResult:
    factor: RedistFactor
    pts: list
    M_I_max: float
    M_II_max: float
    M_inf_max: float
    x_inf_max: float
    M_I_min: float
    M_II_min: float
    M_inf_min: float
    x_inf_min: float

    def at(self, x: float) -> RedistPoint:
        """最接近里程 x 的取樣點。"""
        return min(self.pts, key=lambda p: abs(p.x - x))

    def envelope(self, x: float):
        """該點應設計的彎矩範圍 (min, max)＝合龍當下 M_I 與長期 M(∞) 兩狀態。

        M(t) = M_I + λ(t)(M_II − M_I)，λ(t) 自 0 單調增至 λ∞ < 1 → 極值只會出現在
        兩端點；M_II（一次成形）**實際從未達到**，不必納入。
        逐跨施工的橋既不是簡支也不是連續：跨中保留較多簡支的正彎矩、
        墩頂又拿到大部分連續的負彎矩，兩者都必須設計。
        """
        p = self.at(x)
        return (min(p.M_I, p.M_inf), max(p.M_I, p.M_inf))


def creep_redistribution(xs: Sequence[float], M_I: Sequence[float],
                         M_II: Sequence[float], dphi: float,
                         chi: float = 0.8, method: str = "trost") -> StagingResult:
    """對已知的 M_I／M_II 彎矩圖做潛變重分配。

    xs/M_I/M_II 同長度；M_I 為施工體系解、M_II 為最終體系解（同一組載重）。
    """
    if not (len(xs) == len(M_I) == len(M_II)):
        raise ValueError("xs / M_I / M_II 長度須相同")
    f = redistribution_factor(dphi, chi, method)
    pts = []
    for x, m1, m2 in zip(xs, M_I, M_II):
        d = f.lam * (m2 - m1)
        pts.append(RedistPoint(x=x, M_I=m1, M_II=m2, dM=d, M_inf=m1 + d))
    hi = max(pts, key=lambda p: p.M_inf)
    lo = min(pts, key=lambda p: p.M_inf)
    return StagingResult(
        factor=f, pts=pts,
        M_I_max=max(M_I), M_II_max=max(M_II), M_inf_max=hi.M_inf, x_inf_max=hi.x,
        M_I_min=min(M_I), M_II_min=min(M_II), M_inf_min=lo.M_inf, x_inf_min=lo.x)


# ── 常見工序：逐跨施工（先簡支後連續）──────────────────────────
def simple_span_dl_moment(spans: Sequence[float], w: float, x: float) -> float:
    """各跨獨立簡支時，里程 x 的恆載彎矩 M_I = w·x'(L−x')/2。"""
    acc = 0.0
    for L in spans:
        if x <= acc + L + 1e-9:
            xp = min(max(x - acc, 0.0), L)
            return w * xp * (L - xp) / 2.0
        acc += L
    return 0.0


def span_by_span_dead_load(spans: Sequence[float], w: float, n_per_span: int = 20):
    """逐跨施工（先簡支後連續）的 M_I 與 M_II。

    M_I：各跨簡支（施工時各跨自承）；M_II：同一 w 作用於連續體系。
    回傳 (xs, M_I, M_II)，可直接餵 `creep_redistribution`。
    """
    from .influence_cont import cont_dl_moment
    total = sum(spans)
    n = n_per_span * len(spans)
    xs, m1, m2 = [], [], []
    brk = []
    acc = 0.0
    for L in spans:
        acc += L
        brk.append(acc)
    for i in range(n + 1):
        x = total * i / n
        xs.append(x)
        m1.append(simple_span_dl_moment(spans, w, x))
        m2.append(cont_dl_moment(spans, w, x))
    # 支承處補點：束制彎矩的折點都在支承，漏掉會讓線性檢查與極值失真
    for b in brk[:-1]:
        if all(abs(b - v) > 1e-9 for v in xs):
            xs.append(b)
            m1.append(simple_span_dl_moment(spans, w, b))
            m2.append(cont_dl_moment(spans, w, b))
    order = sorted(range(len(xs)), key=lambda k: xs[k])
    return ([xs[k] for k in order], [m1[k] for k in order], [m2[k] for k in order])


def redistribution_is_linear(res: StagingResult, spans: Sequence[float],
                             tol: float = 1e-6) -> bool:
    """自我驗證：束制彎矩 ΔM 在**支承間必為線性**。

    ΔM = λ(M_II − M_I) 不對應任何外載（兩者是同一組載重的兩個平衡解，差值是
    自平衡的束制場），故其二階微分為零 → 支承間線性。這是本模組最強的檢查：
    只要 M_I 或 M_II 有一個算錯，或取樣點沒對齊支承，這條就會掛掉。
    """
    edges, acc = [0.0], 0.0
    for L in spans:
        acc += L
        edges.append(acc)
    for a, b in zip(edges[:-1], edges[1:]):
        seg = [p for p in res.pts if a - 1e-9 <= p.x <= b + 1e-9]
        if len(seg) < 3:
            continue
        x0, y0 = seg[0].x, seg[0].dM
        x1, y1 = seg[-1].x, seg[-1].dM
        span = x1 - x0
        if span <= 0:
            continue
        scale = max(abs(p.dM) for p in res.pts) or 1.0
        for p in seg:
            lin = y0 + (y1 - y0) * (p.x - x0) / span
            if abs(p.dM - lin) > tol * scale:
                return False
    return True


# ── 預力次彎矩 M₂ 的重分配 ─────────────────────────────────────
@dataclass
class M2RedistResult:
    factor: RedistFactor
    xs: list
    M2_cont: list       # 若一開始就在連續體系張拉，會有的 M₂
    M2_inf: list        # 先簡支張拉後連續，潛變生成的 M₂(∞)
    M2_pier_cont: float
    M2_pier_inf: float


def prestress_M2_redistribution(spans: Sequence[float], xs: Sequence[float],
                                M2_cont: Sequence[float], dphi: float,
                                chi: float = 0.8, method: str = "trost") -> M2RedistResult:
    """先簡支張拉、後建立連續時的次彎矩 M₂(∞)。

    🔴 **這件事很容易被漏掉**：在簡支體系張拉時結構是**靜定**的，次彎矩恆為 0——
      手算到此常就結案。但建立連續後，預力造成的曲率同樣會被新接頭束制，
      潛變會**逐漸生成** M₂，終值為 λ·M₂,cont。

    亦即 M_I = 0、M_II = M₂,cont，套同一個重分配式。

    ⚠ M₂ 與恆載重分配的方向常相反（M₂ 在墩頂多為正、恆載為負），
      故**不可只算其中一個**——兩者要各自重分配後再疊加。
    """
    f = redistribution_factor(dphi, chi, method)
    inf = [f.lam * m for m in M2_cont]
    ip = min(range(len(xs)), key=lambda i: abs(xs[i] - spans[0]))
    return M2RedistResult(factor=f, xs=list(xs), M2_cont=list(M2_cont), M2_inf=inf,
                          M2_pier_cont=M2_cont[ip], M2_pier_inf=inf[ip])


# ── 時序敏感度：連續性建立得多早 ──────────────────────────────
@dataclass
class TimingRow:
    t1_label: str
    phi_t1: float
    dphi: float
    lam: float
    M_pier: float


def timing_sensitivity(phi_inf: float, phi_at_continuity: Sequence,
                       M_I_pier: float, M_II_pier: float,
                       chi: float = 0.8, method: str = "trost") -> list:
    """連續性建立時機對墩頂彎矩的影響。

    phi_at_continuity：[(標籤, φ(t₁,t₀)), ...]——連續時已發生的潛變量。

    🔴 **與直覺相反**：決定重分配的是**剩餘**潛變 Δφ = φ(∞) − φ(t₁)，
      所以**越早建立連續、墩頂負彎矩越大**。趕工早合龍不會讓結構比較輕鬆，
      反而使墩頂更不利；反之延後合龍可降低墩頂負彎矩（代價是跨中正彎矩較高）。
    """
    rows = []
    for label, p1 in phi_at_continuity:
        d = max(0.0, phi_inf - p1)
        f = redistribution_factor(d, chi, method)
        rows.append(TimingRow(t1_label=label, phi_t1=p1, dphi=d, lam=f.lam,
                              M_pier=M_I_pier + f.lam * (M_II_pier - M_I_pier)))
    return rows


# ── 簡支時張拉的腱：各跨自己的拋物線 ────────────────────────────
def simple_span_tendon_segs(spans: Sequence[float], e_mid: float, e_end: float = 0.0):
    """先簡支後連續、**簡支時張拉**的腱線形：各跨獨立拋物線（端 e_end、跨中 e_mid）。

    🔴 不可拿連續腱線形（通過墩頂、有負偏心）來算 M₂,cont——那種線形只有在
      連續**之後**張拉才存在。兩者的 M₂ 差很多（40+40：簡支線形 P·a = 26,311
      vs 連續線形 14,185），2026-09-18 以前的 H7 §六即犯此錯。
    """
    from .continuous import parabola_seg
    segs, acc = [], 0.0
    for L in spans:
        xm = acc + L / 2.0
        segs.append(parabola_seg(acc, acc + L, xm, e_mid, -4.0 * (e_mid - e_end) / (L * L), "simple"))
        acc += L
    return segs


# ── 把體系轉換套進連續梁包絡 ────────────────────────────────────
@dataclass
class StagedEnvRow:
    x: float
    M_dc_I: float       # 合龍當下（施工體系）自重彎矩
    M_dc_II: float      # 一次成形（連續體系）自重彎矩——僅供對照
    M_dc_inf: float     # 長期重分配後
    M_dw: float         # 附加恆載（連續後才施加，不重分配）
    M_ll_pos: float
    M_ll_neg: float
    M2_t1: float        # 合龍當下的預力次彎矩
    M2_inf: float       # 長期
    Ms_pos: float       # 兩狀態取不利
    Ms_neg: float
    Mu_pos: float
    Mu_neg: float
    gov_pos: str        # "t1"／"inf"：正彎矩由哪個狀態控制
    gov_neg: str


def staged_envelope(spans: Sequence[float], rows, w_dc: float, lam: float,
                    M2=None, ps_at_simple: bool = False,
                    g_dc: float = 1.25, g_dw: float = 1.50, g_ll: float = 1.75,
                    g_dc_min: float = 0.90, g_dw_min: float = 0.65):
    """把逐跨施工的潛變重分配套進連續梁包絡（rows＝taiwan_cont_envelope 的列）。

    檢核兩個真實狀態，逐列取不利：
      t₁（合龍當下）：自重 = M_I（各跨簡支）
      ∞ （長期）    ：自重 = M_I + λ(M_II − M_I)
    M(t) 在兩者之間單調變化，極值只在端點；**M_II 從未真的出現**。

    附加恆載 DW 與活載皆在連續後施加 → 直接作用於連續體系、不重分配。

    預力次彎矩 M2（逐列，連續體系下的 M₂,cont）：
      ps_at_simple=False（**連續後才張拉**，如全長連續腱／墩頂頂板腱）：兩狀態皆為全量。
      ps_at_simple=True （**簡支時張拉**）：t₁ 為 0（靜定）、∞ 為 λ·M₂,cont；
                          此時 M2 必須以 `simple_span_tendon_segs` 的線形算。
    M₂ 載重因數 1.0（與分析器一致）。

    ⚠ **載重因數依正負號取 max／min**（AASHTO γ_p）：恆載對所檢核的彎矩**有利**時取
      γ_min（DC 0.90、DW 0.65）。這在墩頂**正彎矩**特別關鍵——自重是負彎矩、對正彎矩
      有利，若仍乘 1.25 會嚴重低估正彎矩接頭的需求。
    AASHTO 5.14.1.4.2「束制彎矩有利時不得計入任何組合」由 t₁（無束制）與 ∞（有束制）
    兩狀態取不利**自動滿足**。
    """
    from .influence_cont import factored   # 與基礎包絡同一套因數規則（避免兩份漂移）
    out = []
    for i, r in enumerate(rows):
        mI = simple_span_dl_moment(spans, w_dc, r.x)
        mII = r.M_dc
        minf = mI + lam * (mII - mI)
        m2 = 0.0 if M2 is None else M2[i]
        m2_t1, m2_inf = (0.0, lam * m2) if ps_at_simple else (m2, m2)
        st = {}
        for tag, dc, mm2 in (("t1", mI, m2_t1), ("inf", minf, m2_inf)):
            st[tag] = (dc + r.M_dw + r.M_ll_pos + mm2, dc + r.M_dw + r.M_ll_neg + mm2,
                       factored(dc, r.M_dw, r.M_ll_pos, True, g_dc, g_dc_min, g_dw, g_dw_min, g_ll) + mm2,
                       factored(dc, r.M_dw, r.M_ll_neg, False, g_dc, g_dc_min, g_dw, g_dw_min, g_ll) + mm2)
        gp = "t1" if st["t1"][2] >= st["inf"][2] else "inf"
        gn = "t1" if st["t1"][3] <= st["inf"][3] else "inf"
        out.append(StagedEnvRow(
            x=r.x, M_dc_I=mI, M_dc_II=mII, M_dc_inf=minf, M_dw=r.M_dw,
            M_ll_pos=r.M_ll_pos, M_ll_neg=r.M_ll_neg, M2_t1=m2_t1, M2_inf=m2_inf,
            Ms_pos=max(st["t1"][0], st["inf"][0]), Ms_neg=min(st["t1"][1], st["inf"][1]),
            Mu_pos=max(st["t1"][2], st["inf"][2]), Mu_neg=min(st["t1"][3], st["inf"][3]),
            gov_pos=gp, gov_neg=gn))
    return out
