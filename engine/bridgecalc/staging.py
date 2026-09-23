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
                          method: str = "trost", phi_restraint: float = None) -> RedistFactor:
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

    phi_restraint（嚴格 AAEM）：相容條件的分子與分母其實是**兩個不同的潛變係數**——
        分子：自重在 t0 加載、t1 之後剩餘的潛變 Δφ = φ(∞,t0) − φ(t1,t0)
        分母：束制彎矩自 t1 起逐漸施加，其潛變為 φ(∞,t1)
        λ = Δφ / (1 + χ·φ(∞,t1))
      不給時以 Δφ 代 φ(∞,t1)（即上式常用的單一 Δφ 近似）。一般 φ(∞,t1) ≥ Δφ，故近似的 λ
      **偏大**——墩頂負彎矩與正束制彎矩偏保守（跨中由 t₁ 狀態控制，不受影響）。
    """
    if dphi < 0:
        raise ValueError("Δφ 不可為負（連續後尚未發生的潛變量）")
    if method == "dischinger":
        lam = 1 - math.exp(-dphi)
    elif method == "trost":
        pr = dphi if phi_restraint is None else phi_restraint
        lam = dphi / (1 + chi * pr) if (1 + chi * pr) != 0 else 0.0
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
                         chi: float = 0.8, method: str = "trost",
                         phi_restraint: float = None) -> StagingResult:
    """對已知的 M_I／M_II 彎矩圖做潛變重分配。

    xs/M_I/M_II 同長度；M_I 為施工體系解、M_II 為最終體系解（同一組載重）。
    """
    if not (len(xs) == len(M_I) == len(M_II)):
        raise ValueError("xs / M_I / M_II 長度須相同")
    f = redistribution_factor(dphi, chi, method, phi_restraint)
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
                                chi: float = 0.8, method: str = "trost",
                                phi_restraint: float = None) -> M2RedistResult:
    """先簡支張拉、後建立連續時的次彎矩 M₂(∞)。

    🔴 **這件事很容易被漏掉**：在簡支體系張拉時結構是**靜定**的，次彎矩恆為 0——
      手算到此常就結案。但建立連續後，預力造成的曲率同樣會被新接頭束制，
      潛變會**逐漸生成** M₂，終值為 λ·M₂,cont。

    亦即 M_I = 0、M_II = M₂,cont，套同一個重分配式。

    ⚠ M₂ 與恆載重分配的方向常相反（M₂ 在墩頂多為正、恆載為負），
      故**不可只算其中一個**——兩者要各自重分配後再疊加。
    """
    f = redistribution_factor(dphi, chi, method, phi_restraint)
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


# ── 體系轉換的恆載剪力 ─────────────────────────────────────────
def simple_span_dl_shear(spans: Sequence[float], w: float, x: float, side: str = "R") -> float:
    """各跨獨立簡支時，斷面 x（side 側）的恆載剪力 V_I = w(L/2 − a)。"""
    from .influence_cont import _supports, _shear_span
    xs = _supports(spans)
    i = _shear_span(xs, spans, x, side)
    return w * (spans[i] / 2.0 - (x - xs[i]))


def staged_shear_row(spans: Sequence[float], row, w_dc: float, lam: float):
    """把一列連續梁剪力（ContShearRow）改為 t₁／∞ 兩狀態取不利。

    t₁：自重剪力＝簡支 V_I；∞：V_I + λ(V_II − V_I)。束制剪力 λ(V_II − V_I) 為束制彎矩之
    斜率，**每跨內為常數**。端支承 t₁ 為 wL/2，大於連續的 3wL/8——只用連續體系會低估。
    回傳同欄位的 ContShearRow（V_dc 為控制狀態值），另附 V_dc_I／V_dc_inf／gov。
    """
    from .influence_cont import ContShearRow, factored
    vI = simple_span_dl_shear(spans, w_dc, row.x, row.side)
    vinf = vI + lam * (row.V_dc - vI)
    cand = []
    for tag, dc in (("t1", vI), ("inf", vinf)):
        cand.append((tag, dc, factored(dc, row.V_dw, row.V_ll_pos, True),
                     factored(dc, row.V_dw, row.V_ll_neg, False)))
    p = max(cand, key=lambda c: c[2])
    n = min(cand, key=lambda c: c[3])
    ctrl = p if abs(p[2]) >= abs(n[3]) else n
    out = ContShearRow(row.x, row.side, ctrl[1], row.V_dw, row.V_ll_pos, row.V_ll_neg,
                       p[2], n[3], row.I)
    out.V_dc_I, out.V_dc_II, out.V_dc_inf = vI, row.V_dc, vinf
    out.gov_pos, out.gov_neg = p[0], n[0]
    return out


# ── 連續橫隔梁正彎矩接頭（AASHTO 5.14.1.4.9，NLM 7839c56d 核對）───────
@dataclass
class PosMomentConnResult:
    fr: float           # 破裂模數 MPa（橫隔梁混凝土）
    Mcr: float          # kN·m
    Mu_pos: float       # 係數化正束制彎矩 kN·m（≤0 表示無）
    simplified: bool    # 合龍時梁齡 ≥ 90 天 → 5.14.1.4.4 簡化
    M_req: float        # 接頭需求彎矩 kN·m
    governs: str        # "Mu+"／"0.6Mcr"／"1.2Mcr"
    As_est: float       # mm²，估算（φ=0.9、jd=0.9d）；d 未給則為 None


def positive_moment_connection(Mu_pos: float, I_g: float, y_t: float, fc_diaph: float,
                               age_days: float = None, fr: float = None,
                               d: float = None, fy: float = 420.0, phi: float = 0.9):
    """連續橫隔梁正彎矩接頭需求（簡支預鑄梁連續化／逐跨施工先簡支後連續）。

    - 5.14.1.4.5：**所有**連續橫隔梁皆須設正、負彎矩接頭，不論連續程度。
    - 5.14.1.4.9a：鋼筋取 max(係數化正束制彎矩, 0.6 M_cr)。
    - 5.14.1.4.4：合龍時梁齡 ≥ 90 天（業主同意、契約規定）→ 正束制可取 0，接頭 ≥ 1.2 M_cr。
    - M_cr＝f_r·I_g/y_t（Eq. 5.7.3.6.2-2）：**總毛複合斷面**、**橫隔梁混凝土**的破裂模數、
      **不計預力**（C5.14.1.4.9a：橫隔梁非預力斷面）。f_r 預設 0.63√f'c（5.4.2.6 基本值；
      0.97√f'c 是最小鋼筋專用，不適用此處）。
    I_g mm⁴、y_t mm（中性軸至受拉緣＝底緣）、fc MPa；As_est 以 M_req/(φ·f_y·0.9d) 估算。
    """
    f_r = 0.63 * math.sqrt(fc_diaph) if fr is None else fr
    Mcr = f_r * I_g / y_t / 1e6
    simp = age_days is not None and age_days >= 90
    if simp:
        M_req, gov = 1.2 * Mcr, "1.2Mcr"
    else:
        M_req, gov = (Mu_pos, "Mu+") if Mu_pos > 0.6 * Mcr else (0.6 * Mcr, "0.6Mcr")
    As = M_req * 1e6 / (phi * fy * 0.9 * d) if d else None
    return PosMomentConnResult(fr=f_r, Mcr=Mcr, Mu_pos=Mu_pos, simplified=simp,
                               M_req=M_req, governs=gov, As_est=As)


# ── 潛變係數 φ(t, t_i)：AASHTO LRFD 5.4.2.3.2（NLM 7839c56d 核對，p.5-15～5-16）──────
# 台灣公路橋梁設計規範**沒有** φ(t) 模式（潛變只以 §8.16.2 經驗式 12f_cir − 7f_cds 處理，
# NLM 7d947294 確認），故體系轉換需要 φ 時採 AASHTO。
@dataclass
class CreepAASHTO:
    ks: float
    khc: float
    kf: float
    ktd: float
    psi: float


def aashto_creep(t: float, ti: float, H: float = 75.0, VS: float = 150.0,
                 fci: float = 32.0) -> CreepAASHTO:
    """AASHTO LRFD Eq. 5.4.2.3.2-1～5（SI）。

        ψ(t,t_i) = 1.9·k_s·k_hc·k_f·k_td·t_i^−0.118
        k_s = 1.45 − 0.0051(V/S) ≥ 1.0    （V/S mm；C5.4.2.3.2：建式時最大 V/S 150 mm）
        k_hc = 1.56 − 0.008H               （H 年平均相對濕度 %）
        k_f = 35/(7 + f'ci)                （f'ci MPa；未知時可取 0.80 f'c）
        k_td = t/(61 − 0.58 f'ci + t)       （t＝**載重後經過**天數；t=∞ → 1）

    t_i：載重時材齡（天）。適用 f'c ≤ 105 MPa（5.4.2.3.1）。
    """
    ks = max(1.0, 1.45 - 0.0051 * VS)
    khc = 1.56 - 0.008 * H
    kf = 35.0 / (7.0 + fci)
    ktd = 1.0 if t == math.inf else t / (61.0 - 0.58 * fci + t)
    return CreepAASHTO(ks=ks, khc=khc, kf=kf, ktd=ktd,
                       psi=1.9 * ks * khc * kf * ktd * ti ** -0.118)


@dataclass
class StagingPhi:
    t0: float           # 載重（自重作用／脫模）材齡 天
    t1: float           # 建立連續材齡 天
    phi_inf: float      # ψ(∞, t0)
    phi_t1: float       # ψ(t1 − t0, t0)：連續前已發生
    dphi_sub: float     # 做法②：全量相減 ψ(∞,t0) − ψ(t1−t0,t0)（Mattock 式）
    phi_load_t1: float  # 做法①：以 t1 為載重材齡 ψ(∞, t1)（便覽 Dischinger 式）


def staging_phi(t0: float, t1: float, H: float = 75.0, VS: float = 150.0,
                fci: float = 32.0) -> StagingPhi:
    """體系轉換用的剩餘潛變 Δφ（AASHTO ψ）。兩種取法並列、不替使用者選邊（見 H7 §七）：

    ② 全量相減：自重在 t0 作用於簡支體系，連續時已發生 ψ(t1−t0, t0)，剩餘＝ψ(∞,t0)−ψ(t1−t0,t0)。
    ① 以轉換材齡查：ψ(∞, t1)——便覽 §17.2.3 以「90 日轉連續取 φ=1.7」即此法。
    兩者差在老化效應；② 對「早加載、晚連續」較貼近實際載重歷程。
    """
    if t1 < t0:
        raise ValueError("連續材齡 t1 不可早於載重材齡 t0")
    inf0 = aashto_creep(math.inf, t0, H, VS, fci).psi
    at1 = aashto_creep(t1 - t0, t0, H, VS, fci).psi if t1 > t0 else 0.0
    return StagingPhi(t0=t0, t1=t1, phi_inf=inf0, phi_t1=at1, dphi_sub=inf0 - at1,
                      phi_load_t1=aashto_creep(math.inf, t1, H, VS, fci).psi)



# ── 箱梁體積表面積比（AASHTO C5.4.2.3.2）──────────────────────────
def box_volume_surface(top_w: float, top_t: float, bot_w: float, bot_t: float,
                       web_t: float, n_web: int, h: float, inner_factor: float = 0.5) -> float:
    """V/S（mm）。C5.4.2.3.2：表面積只計暴露於大氣乾燥者；「For poorly ventilated enclosed
    cells, only 50 percent of the interior perimeter should be used」。
    與 section_from_dims 同一三矩形理想化：腹板垂直、立於底板兩緣。
    """
    hw = h - top_t - bot_t
    A = top_w * top_t + bot_w * bot_t + n_web * web_t * hw
    p_out = top_w + (top_w - bot_w) + 2 * top_t + 2 * hw + 2 * bot_t + bot_w
    p_in = 2 * (bot_w - n_web * web_t) + 2 * max(0, n_web - 1) * hw
    return A / (p_out + inner_factor * p_in)


@dataclass
class TimingRowAASHTO:
    t1: float
    dphi: float          # φ(∞,t0) − φ(t1−t0,t0)
    phi_r: float         # φ(∞,t1)
    lam_exact: float     # Δφ/(1+χ·φ_r)
    lam_approx: float    # Δφ/(1+χ·Δφ)
    M_pier: float        # 以 lam_exact


def timing_sensitivity_aashto(t0: float, t1_list, M_I_pier: float, M_II_pier: float,
                              H: float = 75.0, VS: float = 150.0, fci: float = 32.0,
                              chi: float = 0.8) -> list:
    """以 AASHTO ψ 模式算合龍時機敏感度（取代手填 φ 的示範表）。"""
    out = []
    for t1 in t1_list:
        sp = staging_phi(t0, t1, H, VS, fci)
        ex = redistribution_factor(sp.dphi_sub, chi, "trost", sp.phi_load_t1).lam
        ap = redistribution_factor(sp.dphi_sub, chi, "trost").lam
        out.append(TimingRowAASHTO(t1=t1, dphi=sp.dphi_sub, phi_r=sp.phi_load_t1,
                                   lam_exact=ex, lam_approx=ap,
                                   M_pier=M_I_pier + ex * (M_II_pier - M_I_pier)))
    return out


# ── 正彎矩接頭之錨定與配置（B2，AASHTO 5.14.1.4.9a~d；2026-09-22 NLM 核）──────────
# 5.14.1.4.9a 允許三種接頭：①埋入預鑄梁並伸展進橫隔梁之一般鋼筋；②未於梁端解除握裹之
#   先拉鋼絞線延伸錨入橫隔梁（不得使用 debonded／shielded 者）；③經分析、試驗或業主認可者。
# 5.14.1.4.9b（一般鋼筋）原文重點：
#   - 錨定須符合 Art. 5.11，**且鋼筋須伸展至超過支承面內側邊緣**（beyond the inside edge
#     of the bearing area）——臨界斷面在梁內，不是在橫隔梁面。
#   - 正彎矩鋼筋若加在先拉鋼絞線之間，須考量**混凝土搗實與握裹**。
#   - 多支鋼筋時，**截斷點須成對錯開且對稱於預鑄梁中心線**。
# 5.14.1.4.9c（延伸鋼絞線）：以 90° 彎鉤或依 Art. 5.11.4 之伸展長度錨入橫隔梁；
#   **彎折前自梁面外伸 ≥ 200 mm**；設計應力上限為延伸總長之函數：
#     f_psl = (ℓ_dsh − 203)/0.840 （服務，假設斷面開裂）   Eq. 5.14.1.4.9c-1
#     f_pul = (ℓ_dsh − 203)/0.600 （強度）                Eq. 5.14.1.4.9c-2
#   🔴 式中**沒有**鋼絞線直徑 d_b（英制對應 (ℓ−8)/0.228、(ℓ−8)/0.163 ksi-in，
#      換算 1/(0.228×25.4)×6.895 = 1/0.840 ✓）。ℓ_dsh 為延伸鋼絞線之**總長度**。
# 5.14.1.4.9d：配置須對稱（或儘可能對稱）於斷面中心線；須考量製造與吊裝，
#   兩側梁伸出之鋼筋須能交錯而不衝突，並留出橫隔梁錨定筋之置放空間。
STRAND_PROJECT_MIN = 200.0     # mm，彎折前最小外伸
STRAND_L0 = 203.0              # mm，式中扣除常數（＝8 in.）
STRAND_K_SERVICE = 0.840
STRAND_K_STRENGTH = 0.600


def strand_stress_extended(l_dsh: float):
    """延伸鋼絞線之設計應力上限 (f_psl, f_pul)，MPa。l_dsh：延伸鋼絞線總長度 mm。"""
    eff = max(0.0, l_dsh - STRAND_L0)
    return eff / STRAND_K_SERVICE, eff / STRAND_K_STRENGTH


def strand_length_required(f_target: float, limit: str = "strength") -> float:
    """達到目標應力所需之延伸總長度 ℓ_dsh，mm（上式反解）。"""
    k = STRAND_K_STRENGTH if limit == "strength" else STRAND_K_SERVICE
    return f_target * k + STRAND_L0


@dataclass
class PosConnDetailResult:
    kind: str              # "strand"／"rebar"
    M_req: float           # kN·m
    d: float               # mm，接頭鋼筋至受壓緣之有效深度
    n_req: float           # 需求支數（未取整）
    n_use: int             # 採用支數（鋼絞線取偶數以對稱；鋼筋取偶數成對錯開）
    A_each: float          # mm²
    f_design: float        # MPa，設計應力（鋼絞線為 f_pul；鋼筋為 f_y）
    phiMn: float           # kN·m，採用支數之設計彎矩強度
    ok: bool
    # 鋼絞線專用
    l_dsh: float = None            # mm，延伸總長度
    f_psl: float = None            # MPa，服務狀態上限
    project_ok: bool = None        # 彎折前外伸 ≥200 mm
    l_dsh_req: float = None        # mm，若以 f_y 級應力設計所需長度（對照）
    # 一般鋼筋專用
    ld: float = None               # mm，需求伸展長度（Art. 5.11）
    dev_available: float = None    # mm，自臨界斷面可用之伸展長度
    dev_ok: bool = None            # 是否伸展超過支承面內側邊緣
    stagger_note: str = ""


def pos_conn_strand(M_req: float, d: float, l_dsh: float, A_strand: float = 98.7,
                    projection: float = 250.0, phi: float = 0.9,
                    jd_ratio: float = 0.9) -> PosConnDetailResult:
    """延伸鋼絞線式正彎矩接頭（5.14.1.4.9c）。

    M_req[kN·m] 由 positive_moment_connection() 取得；d[mm]；l_dsh[mm]：延伸總長度；
    A_strand[mm²]：單股面積（15.2 mm 低鬆弛為 140；12.7 mm 為 98.7）；
    projection[mm]：彎折前自梁面外伸長度（須 ≥200）。
    ⚠ 強度以 f_pul 計；**不得使用梁端解除握裹（debonded／shielded）之鋼絞線**。
    """
    f_psl, f_pul = strand_stress_extended(l_dsh)
    arm = jd_ratio * d
    cap_each = phi * f_pul * A_strand * arm / 1e6      # kN·m／股
    n_req = M_req / cap_each if cap_each > 0 else float("inf")
    n_use = int(math.ceil(n_req / 2.0) * 2) if n_req != float("inf") else 0   # 偶數→對稱
    return PosConnDetailResult(
        kind="strand", M_req=M_req, d=d, n_req=n_req, n_use=n_use, A_each=A_strand,
        f_design=f_pul, phiMn=n_use * cap_each, ok=n_use * cap_each >= M_req,
        l_dsh=l_dsh, f_psl=f_psl, project_ok=projection >= STRAND_PROJECT_MIN,
        l_dsh_req=strand_length_required(1860 * 0.5),      # 對照：取 0.5f_pu 所需長度
        stagger_note="配置須對稱於斷面中心線（5.14.1.4.9d），並確認兩側梁鋼筋交錯不衝突")


def pos_conn_rebar(M_req: float, d: float, ld: float, dev_available: float,
                   bar_area: float = 387.0, fy: float = 420.0, phi: float = 0.9,
                   jd_ratio: float = 0.9) -> PosConnDetailResult:
    """一般鋼筋式正彎矩接頭（5.14.1.4.9b）。

    ld[mm]：依 Art. 5.11 之伸展長度需求；dev_available[mm]：自**支承面內側邊緣**起算
    可用之伸展長度——條文要求鋼筋須伸展**超過**該邊緣，臨界斷面在預鑄梁內。
    bar_area：單支面積（D22＝387、D25＝507）。
    """
    arm = jd_ratio * d
    cap_each = phi * fy * bar_area * arm / 1e6
    n_req = M_req / cap_each if cap_each > 0 else float("inf")
    n_use = int(math.ceil(n_req / 2.0) * 2) if n_req != float("inf") else 0
    return PosConnDetailResult(
        kind="rebar", M_req=M_req, d=d, n_req=n_req, n_use=n_use, A_each=bar_area,
        f_design=fy, phiMn=n_use * cap_each, ok=n_use * cap_each >= M_req,
        ld=ld, dev_available=dev_available, dev_ok=dev_available >= ld,
        stagger_note="多支時截斷點須成對錯開且對稱於預鑄梁中心線（5.14.1.4.9b）；"
                     "加在先拉鋼絞線之間時須檢討混凝土搗實與握裹")


# ── 負彎矩接頭（B5，AASHTO 5.14.1.4.5/.6/.7/.8＋5.11.1.2.3；2026-09-22 NLM 核）────────
# 5.14.1.4.5：「Both a positive and negative moment connection … are required for all
#   continuity diaphragms, **regardless of the degree of continuity**.」→ 兩種接頭都要設。
# 5.14.1.4.8 原文重點：
#   - 現場澆置複合橋面板中之鋼筋須依**強度極限**之負設計彎矩配置（5.14.1.4.7 同旨；
#     5.7.3 之規定適用於橋面板與接頭鋼筋）。
#   - 「Longitudinal reinforcement … shall be anchored in regions of the slab that are in
#     compression at strength limit states」——**錨定端必須落在強度狀態下橋面板受壓的區域**，
#     並須符合 5.11.1.2.3；「The termination of this reinforcement shall be staggered.」
#   - 橋面板全部縱向鋼筋皆可計入負彎矩接頭。
#   - 梁間跨越橫隔梁之負彎矩接頭須符合 5.11.5（續接）；有複合橋面板時「permitted」，
#     **無複合橋面板時為 required**。
# 5.11.1.2.3：支承處負彎矩拉力鋼筋**至少 1/3** 須延伸超過反曲點，長度 ≥ max(d, 12d_b, 0.0625·淨跨)。
# 5.14.1.4.6（服務）：損失後若**梁頂於內支承附近出現拉應力**，適用 Table 5.9.4.1.2-1 之限值，
#   但**以 f'c 代入原式中的 f'ci**，並以 Service III 組合計算：
#     有握裹鋼筋（鋼筋應力取 0.5f_y ≤ 210 MPa）：0.63√f'c
#     無握裹鋼筋：min(0.25√f'c, 1.38 MPa)
NEG_ANCHOR_SPAN_FACTOR = 0.0625      # 1/16 淨跨
NEG_ONE_THIRD = 1.0 / 3.0


def neg_top_tension_limit(fc: float, bonded: bool = True) -> float:
    """5.14.1.4.6 之梁頂拉應力限值（Table 5.9.4.1.2-1，f'c 代 f'ci），MPa。"""
    return 0.63 * math.sqrt(fc) if bonded else min(0.25 * math.sqrt(fc), 1.38)


@dataclass
class NegMomentConnResult:
    Mu_neg: float          # kN·m（取正值輸入）
    d: float               # mm
    As_req: float          # mm²
    n_req: float
    n_use: int
    As_provided: float     # mm²（實配；未給則以 n_use 回算）
    ok_strength: bool
    # 錨定（5.11.1.2.3）
    embed_req: float       # mm，max(d, 12d_b, 0.0625·淨跨)
    embed_have: float      # mm
    embed_ok: bool
    n_one_third: int       # 須滿足延伸規定之最少支數（≥1/3 總數）
    # 服務（5.14.1.4.6）
    sigma_top: float       # MPa（拉為正；未給則 None）
    sigma_limit: float
    service_ok: bool
    # 構造
    connection_required: bool   # 無複合橋面板 → 梁間接頭為必須（5.14.1.4.8）
    note: str = ""


def negative_moment_connection(Mu_neg: float, d: float, clear_span: float,
                               bar_area: float = 387.0, db: float = 22.2,
                               embed_beyond_PI: float = None,
                               As_provided: float = None,
                               sigma_top: float = None, fc: float = 40.0,
                               bonded: bool = True, composite_deck: bool = True,
                               fy: float = 420.0, phi: float = 0.9,
                               jd_ratio: float = 0.9) -> NegMomentConnResult:
    """連續橫隔梁負彎矩接頭（橋面板縱向鋼筋）。

    Mu_neg[kN·m]：強度極限之負設計彎矩（**須含時間效應之負束制彎矩**，5.14.1.4.1／.4.2）；
    d[mm]：複合斷面有效深度（橋面板鋼筋至受壓緣）；clear_span[mm]：淨跨；
    embed_beyond_PI[mm]：實際延伸超過反曲點之長度；sigma_top[MPa]：Service III 之梁頂拉應力（拉為正）。
    """
    arm = jd_ratio * d
    As_req = Mu_neg * 1e6 / (phi * fy * arm)
    n_req = As_req / bar_area
    n_use = int(math.ceil(n_req))
    As_prov = As_provided if As_provided is not None else n_use * bar_area
    embed_req = max(d, 12.0 * db, NEG_ANCHOR_SPAN_FACTOR * clear_span)
    have = 0.0 if embed_beyond_PI is None else embed_beyond_PI
    lim = neg_top_tension_limit(fc, bonded)
    return NegMomentConnResult(
        Mu_neg=Mu_neg, d=d, As_req=As_req, n_req=n_req, n_use=n_use, As_provided=As_prov,
        ok_strength=As_prov >= As_req - 1e-9,
        embed_req=embed_req, embed_have=have, embed_ok=have >= embed_req - 1e-9,
        n_one_third=int(math.ceil(n_use * NEG_ONE_THIRD)),
        sigma_top=sigma_top,
        sigma_limit=lim,
        service_ok=(sigma_top is None) or (sigma_top <= lim + 1e-9),
        connection_required=not composite_deck,
        note="錨定端須落在強度狀態下橋面板受壓之區域；截斷點須錯開（5.14.1.4.8）。"
             "橋面板全部縱向鋼筋皆可計入。無複合橋面板時，梁間跨越橫隔梁之接頭為必須（依 5.11.5 續接）。")


# ── 多階段逐跨施工（C1）：各跨材齡不同 ────────────────────────────────
# 現行 creep_redistribution 假設**全橋同時轉為連續**（單一 λ）。逐跨施工實際上是
# 一段一段合龍的：第 k 跨合龍時，先前各跨的自重已在各自的體系上作用了一段時間，
# 其剩餘潛變也各不相同。做法是**逐階段疊加**——每一階段各有自己的 λ：
#
#   M(∞) = Σ_i [ M_I,i + λ_i (M_II,i − M_I,i) ]
#   λ_i = [ψ(∞, t0_i) − ψ(t_c,i − t0_i, t0_i)] / [1 + χ·ψ(∞, t_c,i)]      （嚴格 AAEM）
#
#   M_I,i：第 i 階段之載重作用在**當時體系**上的彎矩
#   M_II,i：同一載重若作用在**最終體系**上的彎矩
#   t0_i：該載重之加載材齡；t_c,i：該階段體系轉換（合龍）之材齡
#
# 🔴 **便覽 §3.5.4 但書**：「材齡差極大（**1 年以上**）之構件相結合時，標準上應考慮施工中
#   產生之斷面力及構件間之材齡差」——超過 365 天時，單一 λ 的簡化不適用，本函式回報 `age_gap_warn`。
# ⚠ 本函式做的是**同一斷面的逐階段疊加**；各階段的 M_I／M_II 由呼叫者以該階段之結構系統算出
#   （簡支、部分連續、全連續各自的彎矩）。工具不替使用者決定施工順序與體系。
AGE_GAP_LIMIT_DAYS = 365.0


@dataclass
class StageSpec:
    name: str
    t0: float          # 該階段載重之加載材齡（天）
    t_c: float         # 該階段體系轉換（合龍）之材齡（天）
    M_I: float         # kN·m，載重作用於當時體系
    M_II: float        # kN·m，同一載重作用於最終體系


@dataclass
class MultiStageResult:
    rows: list             # [(name, lam, M_I, M_II, M_final_i)]
    M_total: float         # kN·m，各階段疊加
    M_I_total: float
    M_II_total: float
    lam_equiv: float       # 等效單一 λ（M_total 落在 M_I_total~M_II_total 之間的比例）
    age_gap: float         # 天，最大合龍材齡差
    age_gap_warn: bool     # > 365 天 → 便覽 §3.5.4 但書
    single_lam: float      # 對照：若以「全部同時連續」單一 λ 計算之結果比例
    M_single: float        # 對照：單一 λ 之總彎矩


def multi_stage_redistribution(stages, chi: float = 0.8, H: float = 75.0,
                               VS: float = 150.0, fci: float = 32.0,
                               single_stage_t0: float = None,
                               single_stage_tc: float = None) -> MultiStageResult:
    """多階段逐跨施工之潛變重分配（逐階段疊加，各階段自有 λ）。

    stages：StageSpec 串列。single_stage_t0／tc：對照用之「全部同時連續」材齡（不給則取
    第一階段之 t0 與**最後**一個 t_c——即「等到最後一跨合龍才一次轉換」的簡化）。
    """
    rows, M_tot, MI_tot, MII_tot = [], 0.0, 0.0, 0.0
    for st in stages:
        sp = staging_phi(st.t0, st.t_c, H, VS, fci)
        lam = redistribution_factor(sp.dphi_sub, chi, "trost", sp.phi_load_t1).lam
        Mi = st.M_I + lam * (st.M_II - st.M_I)
        rows.append((st.name, lam, st.M_I, st.M_II, Mi))
        M_tot += Mi
        MI_tot += st.M_I
        MII_tot += st.M_II
    tcs = [s.t_c for s in stages]
    gap = max(tcs) - min(tcs) if tcs else 0.0
    t0s = single_stage_t0 if single_stage_t0 is not None else stages[0].t0
    tcs1 = single_stage_tc if single_stage_tc is not None else max(tcs)
    sp1 = staging_phi(t0s, tcs1, H, VS, fci)
    lam1 = redistribution_factor(sp1.dphi_sub, chi, "trost", sp1.phi_load_t1).lam
    M_single = MI_tot + lam1 * (MII_tot - MI_tot)
    den = MII_tot - MI_tot
    return MultiStageResult(
        rows=rows, M_total=M_tot, M_I_total=MI_tot, M_II_total=MII_tot,
        lam_equiv=(M_tot - MI_tot) / den if abs(den) > 1e-12 else 0.0,
        age_gap=gap, age_gap_warn=gap > AGE_GAP_LIMIT_DAYS,
        single_lam=lam1, M_single=M_single)


def span_by_span_schedule(n_spans: int, days_per_span: float, t0_offset: float = 28.0,
                          first_cast_age: float = 28.0):
    """由「每跨施工天數」產生各跨之 (t0, t_c)：第 k 跨於第 k·days_per_span 天澆置。

    回傳 [(t0_k, t_c_k)]，其中 t_c 以**最後一跨合龍**為全橋連續時點；
    t0_k＝該跨自重加載材齡（預設澆置後 28 天施拉／落架）。
    ⚠ 這是最單純的等速施工假設；實際進度應由呼叫者直接給 StageSpec。
    """
    out = []
    last = (n_spans - 1) * days_per_span
    for k in range(n_spans):
        cast = k * days_per_span
        t0 = first_cast_age if k == 0 else t0_offset
        # 該跨自加載到全橋合龍之間隔 → 換算成該跨材齡
        t_c = t0 + (last - cast)
        out.append((t0, max(t_c, t0)))
    return out


def stage_moments_span_by_span(spans, w: float, x: float):
    """逐跨施工：把各跨自重對斷面 x 的貢獻拆成各階段之 (M_I,k, M_II,k)，kN·m。

    M_I,k：第 k 跨自重作用於**當時體系**——逐跨施工時該跨為簡支，故只在該跨內有彎矩、
           其他斷面為 0（尚未連續，載重傳不過去）。
    M_II,k：同一載重作用於**最終連續體系**（對單跨加載解支承彎矩後回推）。
    💡 自我驗證：Σ M_I,k ＝ 簡支彎矩、Σ M_II,k ＝ 連續彎矩（疊加原理）——測試逐點比對。
    """
    from .influence_cont import cont_support_moments_uniform, _supports
    xs = _supports(spans)
    res = []
    for k, L in enumerate(spans):
        if xs[k] - 1e-9 <= x <= xs[k + 1] + 1e-9:
            xi = x - xs[k]
            m1 = w * xi * (L - xi) / 2.0
        else:
            m1 = 0.0
        wv = [0.0] * len(spans)
        wv[k] = w
        res.append((m1, _cont_moment_at(spans, xs, wv, x,
                                        cont_support_moments_uniform(spans, wv))))
    return res


def _cont_moment_at(spans, xs, wv, x, M0):
    """連續梁支承彎矩 M0、跨內均佈 wv 時，斷面 x 之彎矩（kN·m）：支承彎矩線性內插＋簡支拋物線。"""
    from .influence_cont import _span_of
    i = _span_of(xs, x)
    L, a = spans[i], x - xs[i]
    return M0[i] + (M0[i + 1] - M0[i]) * a / L + wv[i] * a * (L - a) / 2.0


# ── 懸臂工法之體系轉換（C2）：X₁／X₀ 自動組成 ─────────────────────────────
# H3 算的是「懸臂體系的斷面力」、H7 算的是「體系轉換後的重分配」，中間那一步——
# 把懸臂幾何翻成 M_I／M_II 一對彎矩圖——原本要使用者自己手填。本節把它自動組成。
#
# 🔴 **符號對照**：便覽／H3 卡用 X₀、X₁，本模組（H7）用 M_II、M_I，是同一組東西：
#     X₁ ＝ M_I ＝ 合龍前瞬間**懸臂體系**的斷面力
#     X₀ ＝ M_II ＝ 同一組載重**一次完工**於最終連續體系的斷面力
#     X_final ＝ X₁ + (X₀ − X₁)·λ   ←── 與 M(∞) = M_I + λ(M_II − M_I) 同式
#   （便覽 Dischinger 寫成 (1−e^(−φ))，即 method="dischinger"；Trost 解見 §七）
#
# 🔴 **最容易錯的一步：X₀ 與 X₁ 必須是「同一組載重」。**
#   便覽把 X₀ 寫成「一次完工之死載＋預力」，若照字面把**合龍後才施加**的載重
#   （支架段落架、合龍段自重、鋪裝欄杆 SDL、二期恆載）也算進 X₀ 而 X₁ 沒有，
#   則 (X₀ − X₁) 會含這些載重的**全量**，再乘 λ(<1) → 只拿到 λ 倍，**少算 (1−λ) 倍**。
#   這些載重本來就直接作用在最終體系（M_I = M_II），**不重分配、要全額加**。
#   墩頂為負彎矩者（SDL 幾乎必然）→ 少算的是負彎矩，**偏不安全**。
#   `cantilever_X1X0` 同時回報正確值與此誤用值供對照（`M_lumped_wrong`）。
#
# 幾何模型：每座懸臂墩為一個 `CantileverUnit`（墩心里程＋左右臂長），臂端即合龍點。
#   未被任何懸臂覆蓋的區段視為**支架上澆置**，其自重於合龍後落架時直接上最終體系。
#   兩臂不等長時墩頂彎矩圖有階差 ＝ 傳入墩柱的**不平衡彎矩**（H3 §十）。


@dataclass
class CantileverUnit:
    x_pier: float      # 墩心里程 m
    a_left: float      # 左臂長 m（墩心至左合龍點）
    a_right: float     # 右臂長 m

    @property
    def x_tip_left(self) -> float:
        return self.x_pier - self.a_left

    @property
    def x_tip_right(self) -> float:
        return self.x_pier + self.a_right


def cantilever_units(spans: Sequence[float], piers: Sequence[float],
                     closures: Sequence[float]) -> list:
    """由墩位與合龍點里程組出各懸臂單元；臂端 = 該墩左右最近的合龍點。

    piers：出懸臂之支承里程（通常是內支承）；closures：合龍點里程。
    某側無合龍點時該臂長為 0（該側不出懸臂，例如端跨全在支架上）。
    """
    total = sum(spans)
    for c in closures:
        if not (0.0 <= c <= total):
            raise ValueError(f"合龍點 {c} 不在橋長內")
    units = []
    for xp in piers:
        left = [c for c in closures if c < xp - 1e-9]
        right = [c for c in closures if c > xp + 1e-9]
        a_l = xp - max(left) if left else 0.0
        a_r = min(right) - xp if right else 0.0
        units.append(CantileverUnit(x_pier=xp, a_left=a_l, a_right=a_r))
    for u, v in zip(units[:-1], units[1:]):
        if u.x_tip_right > v.x_tip_left + 1e-9:
            raise ValueError(f"懸臂單元重疊：{u.x_tip_right} > {v.x_tip_left}")
    return units


def cantilever_layout(spans: Sequence[float], end_frac: float = 0.5):
    """由跨徑推出「各內支承出平衡懸臂」之 (piers, closures)。

    內支承之間的跨：兩墩對伸 → 合龍點在該跨**中央**（一條合龍縫）。
    端跨（僅一側出懸臂）：合龍點距墩 `end_frac`×跨長，其餘為**支架段**。
    ⚠ 這只是最常見的排法；實際合龍順序與位置由設計者決定，可直接給 `cantilever_units`。
    """
    xs = [0.0]
    for L in spans:
        xs.append(xs[-1] + L)
    piers = xs[1:-1]
    if not piers:
        return [], []
    closures = []
    if len(spans) >= 1:
        closures.append(piers[0] - end_frac * spans[0])
    for i in range(len(piers) - 1):
        closures.append((piers[i] + piers[i + 1]) / 2.0)
    closures.append(piers[-1] + end_frac * spans[-1])
    return piers, sorted(closures)


def cantilever_covered(units) -> list:
    """各懸臂單元覆蓋之里程區間 [(x1, x2)]（＝合龍前已上結構的自重範圍）。"""
    return [(u.x_tip_left, u.x_tip_right) for u in units]


def cantilever_falsework(spans: Sequence[float], units) -> list:
    """未被懸臂覆蓋之區間（支架上澆置，落架時直接作用於最終體系）。"""
    total = sum(spans)
    out, cur = [], 0.0
    for x1, x2 in sorted(cantilever_covered(units)):
        if x1 > cur + 1e-9:
            out.append((cur, x1))
        cur = max(cur, x2)
    if cur < total - 1e-9:
        out.append((cur, total))
    return out


def cantilever_M_I(units, w: float, x: float, side: str = "R") -> float:
    """懸臂體系（靜定）下自重 w 於斷面 x 之彎矩 kN·m（負＝頂緣受拉）。

    自由體取**外伸側**：右臂 M = −w(x_tip_R − x)²/2、左臂 M = −w(x − x_tip_L)²/2。
    ⚠ x 恰在墩心且兩臂不等長時彎矩圖有階差（差額傳入墩柱），由 `side` 指定取哪一側。
    """
    for u in units:
        if u.x_tip_left - 1e-9 <= x <= u.x_tip_right + 1e-9:
            if abs(x - u.x_pier) < 1e-9:
                a = u.a_right if side.upper().startswith("R") else u.a_left
                return -w * a * a / 2.0
            if x > u.x_pier:
                return -w * (u.x_tip_right - x) ** 2 / 2.0
            return -w * (x - u.x_tip_left) ** 2 / 2.0
    return 0.0


def cantilever_M_I_points(units, loads, x: float, side: str = "R",
                          overhang: float = 0.0) -> float:
    """懸臂體系下**集中載重** loads=[(p, P)] 於斷面 x 之彎矩 kN·m。

    只有與 x 同臂、且比 x 更外側的載重才對 x 產生彎矩（內側者走支承傳力）。
    🔴 `overhang`：**掛籃（Form Traveler）永遠掛在已澆節塊尖端之外**，若以單元臂長
      硬性過濾就會被整個丟掉（掛籃項在懸臂初期還大於節塊自重，丟掉是嚴重低估）。
      臂端外 `overhang` 公尺內之載重仍歸該臂承擔——檢核最大懸臂狀態時，臂長取**該
      步驟已澆範圍**，掛籃伸出量由此參數帶入。
    """
    m = 0.0
    for u in units:
        if not (u.x_tip_left - 1e-9 <= x <= u.x_tip_right + 1e-9):
            continue
        right = (x > u.x_pier + 1e-9) or (abs(x - u.x_pier) < 1e-9
                                          and side.upper().startswith("R"))
        lo, hi = u.x_tip_left - overhang, u.x_tip_right + overhang
        for p, P in loads:
            if not (lo - 1e-9 <= p <= hi + 1e-9):
                continue
            if right and p >= x - 1e-9:
                m += -P * (p - x)
            elif (not right) and p <= x + 1e-9:
                m += -P * (x - p)
        return m
    return 0.0


def cantilever_unbalanced(units, w: float, extras=None, overhang: float = 0.0) -> list:
    """各墩之不平衡彎矩 kN·m（＝墩心兩側彎矩階差，傳入墩柱）。

    w：自重；extras=[(p, P)] 供掛籃、不對稱節塊等集中載重。
    兩臂等長且無 extras 時恆為 0——這是本函式最直接的自我驗證。
    """
    extras = list(extras or [])
    out = []
    for u in units:
        mL = cantilever_M_I([u], w, u.x_pier, "L") + cantilever_M_I_points([u], extras, u.x_pier, "L", overhang)
        mR = cantilever_M_I([u], w, u.x_pier, "R") + cantilever_M_I_points([u], extras, u.x_pier, "R", overhang)
        out.append((u.x_pier, mL - mR))
    return out


def _cant_grid(spans, units, n_per_span: int = 20):
    """取樣格點：均分格 ＋ 支承 ＋ 墩心 ＋ 合龍點（折點漏掉會使線性檢查失真）。"""
    total = sum(spans)
    n = n_per_span * len(spans)
    xs = [total * i / n for i in range(n + 1)]
    key = []
    acc = 0.0
    for L in spans:
        acc += L
        key.append(acc)
    for u in units:
        key += [u.x_pier, u.x_tip_left, u.x_tip_right]
    for k in key:
        if 0.0 <= k <= total and all(abs(k - v) > 1e-9 for v in xs):
            xs.append(k)
    return sorted(xs)


def cantilever_dead_load(spans: Sequence[float], units, w: float,
                         n_per_span: int = 20, side: str = "R"):
    """懸臂工法自重之 (xs, M_I, M_II)，可直接餵 `creep_redistribution`。

    M_I：懸臂體系（合龍前瞬間）＝ X₁；只有懸臂覆蓋段的自重已上結構。
    M_II：**同一組載重**（同樣只有懸臂覆蓋段）作用於最終連續體系 ＝ X₀。
    ⚠ 支架段自重與合龍後載重不在此列（見模組註解：混進來會被乘 λ 而少算）。
    """
    from .influence_cont import cont_moment_partial_udl
    loads = [(x1, x2, w) for x1, x2 in cantilever_covered(units)]
    xs = _cant_grid(spans, units, n_per_span)
    m1 = [cantilever_M_I(units, w, x, side) for x in xs]
    m2 = [cont_moment_partial_udl(spans, loads, x) for x in xs]
    return xs, m1, m2


@dataclass
class CantileverXResult:
    lam: float
    xs: list
    M_I: list              # X₁（懸臂體系）
    M_II: list             # X₀（同載重一次完工）
    M_inf: list            # X_final = X₁ + λ(X₀ − X₁)
    piers: list            # [(x, X1, X0, X_final)]
    closures: list         # [(x, X1, X0, X_final)]
    unbalanced: list       # [(x_pier, M_unbal)]
    linear_ok: bool        # 束制彎矩 ΔM 支承間線性（自我驗證）
    M_post: list           # 合龍後載重（支架段落架＋SDL）於最終體系之彎矩
    M_dead_total: list     # 正解：M_inf + M_post
    M_lumped_wrong: list   # 誤用：X₁ + λ(X₀_全載重 − X₁)
    lumped_err_pier: float # 兩者在（第一座）墩頂之差 kN·m


def cantilever_X1X0(spans: Sequence[float], units, w: float, dphi: float,
                    w_sdl: float = 0.0, chi: float = 0.8, method: str = "trost",
                    phi_restraint: float = None, n_per_span: int = 20,
                    side: str = "R") -> CantileverXResult:
    """懸臂工法 X₁／X₀ 自動組成＋潛變重分配（H3 → H7）。

    w：上部結構自重 kN/m；w_sdl：合龍後施加之二期恆載 kN/m（全長）。
    支架段自重（未被懸臂覆蓋者）自動歸入「合龍後載重」——落架時體系已連續。
    """
    from .influence_cont import cont_moment_partial_udl, cont_dl_moment
    xs, m1, m2 = cantilever_dead_load(spans, units, w, n_per_span, side)
    res = creep_redistribution(xs, m1, m2, dphi, chi, method, phi_restraint)
    lam = res.factor.lam
    post_loads = [(x1, x2, w) for x1, x2 in cantilever_falsework(spans, units)]
    m_post = [cont_moment_partial_udl(spans, post_loads, x)
              + (cont_dl_moment(spans, w_sdl, x) if w_sdl else 0.0) for x in xs]
    m_inf = [p.M_inf for p in res.pts]
    m_tot = [a + b for a, b in zip(m_inf, m_post)]
    # 誤用對照：把合龍後載重一併塞進 X₀
    m2_all = [a + b for a, b in zip(m2, m_post)]
    m_wrong = [a + lam * (b - a) for a, b in zip(m1, m2_all)]

    def pick(x):
        i = min(range(len(xs)), key=lambda k: abs(xs[k] - x))
        return (xs[i], m1[i], m2[i], m_inf[i])

    piers = [pick(u.x_pier) for u in units]
    cls = sorted({u.x_tip_left for u in units} | {u.x_tip_right for u in units})
    closures = [pick(c) for c in cls
                if any(abs(c - u.x_pier) > 1e-9 for u in units)]
    ip = min(range(len(xs)), key=lambda k: abs(xs[k] - units[0].x_pier))
    return CantileverXResult(
        lam=lam, xs=xs, M_I=m1, M_II=m2, M_inf=m_inf,
        piers=piers, closures=closures,
        unbalanced=cantilever_unbalanced(units, w),
        linear_ok=redistribution_is_linear(res, spans),
        M_post=m_post, M_dead_total=m_tot, M_lumped_wrong=m_wrong,
        lumped_err_pier=m_wrong[ip] - m_tot[ip])


def cantilever_schedule(n_seg: int, days_per_seg: float, t0: float = 7.0,
                        days_to_closure: float = 0.0):
    """由節塊施工循環產生各節塊之 (t0, t_c)：第 k 塊（k=0 最靠墩）於 k·days 澆置。

    t0：掛籃移位（自重上結構）時之材齡；合龍日 ＝ n·days_per_seg + days_to_closure。
    → 越晚澆的節塊合龍時材齡越小、剩餘潛變越多、λ 越大（與逐跨施工同向）。
    """
    close_day = n_seg * days_per_seg + days_to_closure
    out = []
    for k in range(n_seg):
        t_c = close_day - k * days_per_seg
        out.append((t0, max(t_c, t0)))
    return out


def cantilever_stages(spans: Sequence[float], units, segments, x: float,
                      schedule=None, side: str = "R", overhang: float = 0.0):
    """把各節塊自重組成 `StageSpec` 串列（餵 `multi_stage_redistribution`）。

    segments：[(name, p, G)]＝節塊名、形心里程 m、自重 kN（**同一臂由內而外排序**）。
    schedule：[(t0, t_c)] 與 segments 等長；不給時以 `cantilever_schedule` 之預設。
    M_I,k：該節塊自重於懸臂體系對斷面 x 之彎矩；M_II,k：同一集中載重於最終連續體系。
    💡 自我驗證：Σ(M_II,k − M_I,k) 沿梁在支承間線性（束制場自平衡）。
    """
    from .influence_cont import cont_moment_il
    segs = list(segments)
    sch = list(schedule) if schedule is not None else cantilever_schedule(len(segs), 7.0)
    if len(sch) != len(segs):
        raise ValueError("schedule 長度須與 segments 相同")
    out = []
    for (name, p, G), (t0, tc) in zip(segs, sch):
        mi = cantilever_M_I_points(units, [(p, G)], x, side, overhang)
        mii = G * cont_moment_il(spans, x, p)
        out.append(StageSpec(name=name, t0=t0, t_c=tc, M_I=mi, M_II=mii))
    return out
