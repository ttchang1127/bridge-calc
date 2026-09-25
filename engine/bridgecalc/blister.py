"""中間錨碇齒塊（Blister／Anchor Rib）錨碇區設計。

對應 算例_外置PT補強_錨固齒塊中間錨碇設計、公式卡_STM壓拉桿設計。

★ 齒塊與端部錨碇、與轉向塊都不同，差在**後方邊界沒有鋼腱延伸提供回力**：
   - 端部錨碇：力打進端橫隔版，後方就是梁端自由面，靠 bursting 配筋即可。
   - 轉向塊  ：外力只有偏向力 P·sinα（5° 時約 P 的 9%），鋼腱繼續延伸提供回力。
   - 中間錨碇齒塊：外力是**完整的 P_s**，且後方無回力 → **Tie-back 配筋是唯一
     抵抗後向分離的機制，任何情況下不得省略**；介面剪力摩擦亦逼近全腱力。
   四類配筋（爆裂／Tie-back／剝裂／介面剪力摩擦）缺一不可。

單位：力 kN、長度 mm、應力 MPa、鋼筋面積 mm²。
"""
import math
from dataclasses import dataclass


# ── 局部區（Local Zone）：錨板承壓 ──────────────────────────────
@dataclass
class LocalBearingResult:
    f_b: float          # 錨板承壓應力 MPa
    f_b_allow: float    # 容許承壓 MPa（已取上限）
    f_b_uncapped: float # 未套 1.5f'ci 上限前的值（診斷用）
    capped: bool        # 是否由 1.5f'ci 上限控制
    ratio: float        # f_b / f_b_allow
    ok: bool            # 不配螺旋筋是否就通過
    spiral_factor_req: float  # 需要螺旋筋提供幾倍承壓提升


def blister_local_bearing(Pd_kN: float, A_plate: float, A_bearing: float,
                          fci: float, phi_b: float = 0.70) -> LocalBearingResult:
    """錨板局部承壓驗核（AASHTO 局部區）。

        f_b = P_d / A_p
        f_b,allow = 0.85·φ_b·f'ci·√(A_bearing/A_p) ≤ 1.5·f'ci

    Pd_kN：局部驗核用設計力（張拉最不利工況，常取 1.2·P_j）；
    A_plate：錨板面積 mm²；A_bearing：含錨板四周均布影響面積 mm²；
    f'ci：張拉時混凝土強度 MPa。

    ⚠ 齒塊的錨板面積受齒塊尺寸限制，**常常算不過**——這不是錯誤而是設計訊息：
      必須配錨具螺旋筋（可提升 2～3 倍）或改用承壓板較大的錨具系統。
    """
    f_b = Pd_kN * 1e3 / A_plate
    raw = 0.85 * phi_b * fci * math.sqrt(A_bearing / A_plate)
    cap = 1.5 * fci
    allow = min(raw, cap)
    return LocalBearingResult(
        f_b=f_b, f_b_allow=allow, f_b_uncapped=raw, capped=raw > cap,
        ratio=f_b / allow, ok=f_b <= allow,
        spiral_factor_req=max(1.0, f_b / allow))


# ── (a) 爆裂力配筋（Bursting）──────────────────────────────────
@dataclass
class BurstResult:
    F_burst: float      # kN
    d_burst: float      # mm，合力位置（距錨板面）
    As_burst: float     # mm²
    zone_from: float    # mm，配置範圍起（0.5 d_burst）
    zone_to: float      # mm，配置範圍迄（1.5 d_burst）


def blister_bursting(Pu_kN: float, a_plate: float, h: float, e: float = 0.0,
                     fy: float = 420.0, phi: float = 0.9) -> BurstResult:
    """爆裂力與配筋（AASHTO Art. 5.10.9.6.3 一般區簡化式）。

        F_burst = 0.25·P_u·(1 − a/h)
        d_burst = 0.5·(h − 2e)

    a_plate：錨板在該方向的寬度；h：齒塊該方向的總深度；e：錨板對該斷面形心的偏心。
    """
    F = 0.25 * Pu_kN * (1 - a_plate / h)
    d = 0.5 * (h - 2 * abs(e))
    As = F * 1e3 / (phi * fy)
    return BurstResult(F_burst=F, d_burst=d, As_burst=As,
                       zone_from=0.5 * d, zone_to=1.5 * d)


# ── (b) Tie-back 配筋（向後拉結）— 齒塊的關鍵項 ─────────────────
@dataclass
class TiebackResult:
    T_req: float        # kN，0.25·P_s
    C_precomp: float    # kN，既有預壓可抵扣的部分 f_cb·A_cb
    fs_allow: float     # MPa，min(0.6fy, 248)
    As: float           # mm²，含預壓抵扣
    As_conservative: float  # mm²，不抵扣（補強既有橋時預壓分布難確認）
    max_dist_from_axis: float  # mm，配筋須落在距腱軸此距離內（＝一個錨板寬）


def blister_tieback(Ps_kN: float, f_cb: float = 0.0, A_cb: float = 0.0,
                    fy: float = 420.0, a_plate: float = 200.0) -> TiebackResult:
    """Tie-back 配筋（AASHTO Art. 5.10.9.3.4b）。

        A_tie = (0.25·P_s − f_cb·A_cb) / f_s,allow,　f_s,allow = min(0.6 f_y, 248 MPa)

    P_s：**服務狀態**（非係數化）鋼腱力——規範明文以未係數化力的 25% 為下限；
    f_cb：齒塊後方腹板的平均預壓力 MPa；A_cb：齒塊後方板面積 mm²。

    ⚠ 補強既有橋時 f_cb 往往無法確認（既有預力已損失多少、分布如何），此時應取
      `As_conservative`（不抵扣）。本引擎兩者都回傳，不替使用者選邊。
    """
    T = 0.25 * Ps_kN
    C = f_cb * A_cb / 1e3
    fs = min(0.6 * fy, 248.0)
    return TiebackResult(
        T_req=T, C_precomp=C, fs_allow=fs,
        As=max(0.0, (T - C)) * 1e3 / fs,
        As_conservative=T * 1e3 / fs,
        max_dist_from_axis=a_plate)


# ── (c) 剝裂力配筋（Spalling）─────────────────────────────────
@dataclass
class SpallResult:
    F_spall: float
    As_spall: float


def blister_spalling(Ps_kN: float, fy: float = 420.0, phi: float = 0.9,
                     ratio: float = 0.02) -> SpallResult:
    """角隅剝裂力與配筋（現行：0.02·傳入力，blister_design 傳 P_s）。

    🔴 **2026-09-26 稽核發現（待裁示，未改行為）**：兩制皆規定以**乘因數**力計——
      台灣 §8.21.3 4.(8)（p.170）「剝裂力應不小於**乘因數**預力之2%」；
      AASHTO 5.10.9.3.2（p.5-138）「2 percent of the total **factored** tendon force」
      （錨碇區設計力＝1.2×最大張拉力，3.4.3.2）。原註「Art. 5.10.9.6.5，約 0.02·P_s」出處與基準皆不符。
      `blister_design` 目前傳 P_s，B1 算例少算 41%（P_u/P_s＝1.70）。更正會動 golden blister_B1。
    """
    F = ratio * Ps_kN
    return SpallResult(F_spall=F, As_spall=F * 1e3 / (phi * fy))


# ── (d) 介面剪力摩擦配筋（Interface Shear Friction）────────────
# 介面剪應力上限 min(K₁·f'c, K₂) 的兩制常數（2026-09-23 稽核線 NLM 核對）：
#   AASHTO 5.8.4.1（決策 25）：K₁=0.25、K₂=10.3 MPa（整澆與鑿毛面；現澆橋面板於鑿毛預鑄梁頂為 0.30/12.4）
#   台灣公路橋梁設計規範 §7.3.6 4.(4)d（p.135）：「V_n 不得取大於 0.2f'c·A_cv 或 56.2A_cv (kgf/cm²)
#     (5.52A_cv (MPa))」→ 等效 K₁=0.20、K₂=5.52 MPa
# 🔴 **台灣上限遠嚴於 AASHTO**：f'c=40 MPa 時 AASHTO min(10.0,10.3)=10.0、台灣 min(8.0,5.52)=**5.52**
#   （僅 55%）。齒塊／節塊介面以 AASHTO 值通過者，依台灣規範**未必通過**。
# ⚠ 台灣μ 四值與建築規範 112 §22.9 相同（1.4λ/1.0λ/0.6λ/0.7λ），但**上限式不同**——
#   建築分 A/B 兩級（A 級另有 33.6+0.08f'c 與 112 kgf/cm²），橋梁只有單一上限，
#   f'c≥350 kgf/cm² 起建築式偏高（f'c=560 時 +39%），橋梁案不可套用。
K_INTERFACE_AASHTO = (0.25, 10.3)
K_INTERFACE_TW = (0.20, 5.52)

@dataclass
class InterfaceResult:
    V_int: float        # kN，介面剪力
    As_vf: float        # mm²
    tau: float          # MPa，介面平均剪應力
    tau_cap: float      # MPa，介面剪應力上限 min(K1·f'c, K2)
    V_cap: float        # kN，該介面面積下的剪力上限
    area_ok: bool       # 介面面積是否足夠（配筋無法救的那一項）
    A_req: float        # mm²，達到上限所需的最小介面面積


def blister_interface_shear(Ps_kN: float, alpha_deg: float = 5.0,
                            mu: float = 1.0, fsd: float = 360.0,
                            Nd_kN: float = 0.0, gamma0: float = 1.0,
                            A_interface: float = 0.0, fc: float = 40.0,
                            K1: float = 0.25, K2: float = 10.3) -> InterfaceResult:
    """齒塊與腹板介面的剪力摩擦配筋與**面積上限**檢核。

        V_d = P_s·cos α,　γ₀·V_d ≤ μ·(f_sd·A_vf − N_d)
        且  V_ni ≤ min(K₁·f'c, K₂)·A_cv　（配筋再多也不能超過的上限）

    α：鋼腱在錨碇處與橋軸的偏折角（度）。μ：介面粗糙度係數（鑿毛 1.0）。

    ⚠ 齒塊的介面剪力**逼近全腱力**（α 小時 cos α≈1），與轉向塊只傳 P·sinα 差約一個數量級——
      這是齒塊配筋量遠大於轉向塊的主因。

    🔴 **上限檢核不可略**：剪力摩擦的配筋量與介面面積是兩件事——`As_vf` 算得出來不代表
      做得到。介面面積不足時**加多少鋼筋都沒用**，只能加大齒塊。算例常只算配筋而漏掉這項。

    K1/K2 隨規範不同：預設 `K_INTERFACE_AASHTO`＝(0.25, 10.3)；台灣公路橋梁規範 §7.3.6 4.(4)d
    為 `K_INTERFACE_TW`＝(0.20, 5.52)，**遠嚴於 AASHTO**（f'c=40 時 5.52 vs 10.0，僅 55%）。
    採台灣規範時請明確傳入 `K1, K2 = K_INTERFACE_TW`。
    """
    V = Ps_kN * math.cos(math.radians(alpha_deg))
    As = max(0.0, gamma0 * V * 1e3 / mu + Nd_kN * 1e3) / fsd
    tau_cap = min(K1 * fc, K2)
    V_cap = tau_cap * A_interface / 1e3 if A_interface else 0.0
    return InterfaceResult(
        V_int=V, As_vf=As,
        tau=(V * 1e3 / A_interface) if A_interface else 0.0,
        tau_cap=tau_cap, V_cap=V_cap,
        area_ok=(A_interface > 0 and V <= V_cap + 1e-9),
        A_req=V * 1e3 / tau_cap if tau_cap > 0 else float("inf"))


# ── 幾何檢核：錨具承壓面起的直線段與錨板邊距 ─────────────────
@dataclass
class BlisterGeomResult:
    straight_req: float     # mm，需求直線段長
    straight_have: float    # mm
    straight_ok: bool
    edge_req: float         # mm，錨板至齒塊邊的最小邊距
    edge_have: float        # mm
    edge_ok: bool
    L_min: float            # mm，齒塊沿橋軸最小長度
    L_ok: bool
    ok: bool


def blister_geometry_check(L_b: float, W_b: float, a_plate: float, b_plate: float,
                           straight_have: float, straight_req: float = 400.0,
                           edge_req: float = 50.0) -> BlisterGeomResult:
    """齒塊尺寸檢核。

    straight_req：錨具承壓面起必須維持的**直線段**長度（PTI／錨具廠商要求，
    典型 400 mm）——鋼腱不可一離開錨板就起彎，否則錨板受偏心加載且局部承壓惡化。
    edge_req：錨板外緣至齒塊邊的淨距。
    L_b/W_b：齒塊沿橋軸長度／突出腹板面寬度；a_plate×b_plate：錨板尺寸。
    """
    L_min = a_plate + 2 * edge_req
    edge_have = (W_b - b_plate) / 2.0
    s_ok = straight_have >= straight_req - 1e-9
    e_ok = edge_have >= edge_req - 1e-9
    l_ok = L_b >= L_min - 1e-9
    return BlisterGeomResult(
        straight_req=straight_req, straight_have=straight_have, straight_ok=s_ok,
        edge_req=edge_req, edge_have=edge_have, edge_ok=e_ok,
        L_min=L_min, L_ok=l_ok, ok=s_ok and e_ok and l_ok)


# ── 彙整：一次跑完四類配筋＋局部承壓 ──────────────────────────
@dataclass
class BlisterDesign:
    bearing: LocalBearingResult
    burst: BurstResult
    tie: TiebackResult
    spall: SpallResult
    face: InterfaceResult
    geom: BlisterGeomResult
    As_total_conservative: float   # mm²，四類保守值合計（配筋量級參考）
    governing: str                 # 配筋量最大者


def blister_design(Ps_kN: float, Pu_kN: float, Pd_kN: float,
                   a_plate: float = 200.0, b_plate: float = 200.0,
                   L_b: float = 400.0, W_b: float = 300.0, D_b: float = 250.0,
                   A_bearing: float = 60000.0, fci: float = 35.0,
                   f_cb: float = 0.0, A_cb: float = 0.0,
                   alpha_deg: float = 5.0, mu: float = 1.0,
                   fy: float = 420.0, fsd: float = 360.0, fc: float = 40.0,
                   straight_have: float = 400.0) -> BlisterDesign:
    """中間錨碇齒塊完整設計（四類配筋＋局部承壓＋幾何）。

    三個設計力刻意分開傳，因為**它們的來源與係數都不同**，混用是常見錯誤：
      P_s：服務狀態鋼腱力 → Tie-back（規範明文以未係數化力計）與剝裂、介面剪力
      P_u：係數化鋼腱力   → 爆裂（一般區 STM）
      P_d：張拉工況設計力 → 局部區承壓（螺旋筋設計）
    """
    A_plate = a_plate * b_plate
    A_int = L_b * D_b
    bearing = blister_local_bearing(Pd_kN, A_plate, A_bearing, fci)
    burst = blister_bursting(Pu_kN, a_plate, L_b, 0.0, fy)
    tie = blister_tieback(Ps_kN, f_cb, A_cb, fy, a_plate)
    spall = blister_spalling(Ps_kN, fy)
    face = blister_interface_shear(Ps_kN, alpha_deg, mu, fsd, 0.0, 1.0, A_int, fc)
    geom = blister_geometry_check(L_b, W_b, a_plate, b_plate, straight_have)
    items = {"爆裂": burst.As_burst, "Tie-back": tie.As_conservative,
             "剝裂": spall.As_spall, "介面剪力摩擦": face.As_vf}
    return BlisterDesign(
        bearing=bearing, burst=burst, tie=tie, spall=spall, face=face, geom=geom,
        As_total_conservative=sum(items.values()),
        governing=max(items, key=items.get))


# ── 齒塊 STM 壓桿與節點（B1，AASHTO 5.6.3＋5.10.9.3.4c；2026-09-21 NLM 核）─────
# 5.10.9.3.4c 原文：「Reinforcement shall be provided throughout blisters or ribs as required
#   for shear friction, corbel action, bursting forces, and deviation forces due to tendon
#   curvature.」——齒塊本體除既有四類配筋外，**壓桿與節點承載力**須另行驗核（corbel action）。
# ⚠ **STM 的拓樸（壓桿走向、節點型式、拉桿位置）是設計者的選擇**，本函式不替使用者選：
#   壓桿力 C 與夾角 θ 由呼叫者給定；只負責算幾何寬度、有效強度與判定。
# 壓桿寬（CCT 節點，公式卡 F2 §3.2）：w_s = l_b·sinθ + w_t·cosθ
#   l_b＝承壓長度（錨板沿節點面尺寸）、w_t＝拉桿寬度、θ＝壓桿與拉桿夾角。
@dataclass
class BlisterSTMResult:
    C_strut: float        # 壓桿力 N
    theta_deg: float
    w_s: float            # 壓桿有效寬 mm
    b_eff: float          # 壓桿有效厚（取 min(錨板寬, 腹板厚)）mm
    A_cs: float           # 壓桿有效斷面 mm²
    fcu_aashto: float     # AASHTO 2008 式 f_cu MPa
    phiFns_aashto: float  # N
    strut_ok_aashto: bool
    fcu_aci: float        # ACI β 表對照 MPa
    phiFns_aci: float
    strut_ok_aci: bool
    A_n: float            # 節點承壓面積 mm²（預設＝錨板×有效厚；根部節點應另給）
    A_n_req: float        # 達到 φF_nn ≥ P_node 所需之節點面積 mm²（AASHTO 式）
    A_cs_req: float       # 達到 φF_ns ≥ C 所需之壓桿斷面 mm²（AASHTO 式）
    node_type: str
    phiFnn_aashto: float
    node_ok_aashto: bool
    phiFnn_aci: float
    node_ok_aci: bool
    eps1: float           # 主拉應變（AASHTO 式）
    util_strut: float     # 壓桿利用率（AASHTO）
    util_node: float


def blister_stm(C_strut_kN: float, P_node_kN: float, a_plate: float = 200.0,
                b_plate: float = 200.0, w_tie: float = 150.0, theta_deg: float = 45.0,
                web_t: float = 350.0, fc: float = 40.0, node_type: str = "CCT",
                eps_s: float = None, beta_s_aci: str = "reinforced",
                phi: float = 0.70, A_node: float = None) -> BlisterSTMResult:
    """齒塊壓桿與節點驗核；AASHTO 2008 式與 ACI β 表**兩案並列**（決策 19：不替使用者選規範版本）。

    C_strut_kN：壓桿力（係數化）；P_node_kN：節點承壓力（通常＝係數化錨碇力）；
    w_tie：拉桿有效寬（CCT 節點背面深度，＝2×外側面至拉桿重心）；θ：壓桿與拉桿夾角。
    A_node：節點承壓面積，預設＝錨板 × 有效厚。

    🔴 **錨板面屬「局部區」（local zone），不適用一般區之節點限值**——局部區依 AASHTO
    5.10.9.7（承壓試驗／局部承壓＋螺旋圍束）設計，見 `blister_local_bearing`。
    本檢核針對**齒塊根部節點**（力進入腹板處），請以 `A_node` 給實際根部面積（如 L_b × 腹板厚）。
    結果另回傳 `A_n_req`／`A_cs_req`——面積不足時**加鋼筋救不了**，只能加大齒塊；
    AASHTO 5.10.9.3.4a 亦明示「whenever practical … 擴成連續肋（continuous rib）」。
    """
    from .stm import (strut_fcu_aashto, node_capacity_aashto, f_cu, BETA_STRUT,
                      node_capacity, strut_eps1, EPS_S_YIELD, NODE_LIMIT_AASHTO)
    es = EPS_S_YIELD if eps_s is None else eps_s
    th = math.radians(theta_deg)
    w_s = a_plate * math.sin(th) + w_tie * math.cos(th)
    b_eff = min(b_plate, web_t)
    A_cs = w_s * b_eff
    C = C_strut_kN * 1e3
    Pn = P_node_kN * 1e3
    fcu_a = strut_fcu_aashto(fc, es, theta_deg)
    phiFns_a = phi * fcu_a * A_cs
    fcu_i = f_cu(fc, BETA_STRUT[beta_s_aci])
    phiFns_i = phi * fcu_i * A_cs
    A_n = a_plate * b_eff if A_node is None else A_node
    phiFnn_a = node_capacity_aashto(fc, node_type, A_n, phi)
    phiFnn_i = node_capacity(fc, node_type, A_n, phi)
    return BlisterSTMResult(
        C_strut=C, theta_deg=theta_deg, w_s=w_s, b_eff=b_eff, A_cs=A_cs,
        fcu_aashto=fcu_a, phiFns_aashto=phiFns_a, strut_ok_aashto=phiFns_a >= C,
        fcu_aci=fcu_i, phiFns_aci=phiFns_i, strut_ok_aci=phiFns_i >= C,
        A_n=A_n,
        A_n_req=Pn / (NODE_LIMIT_AASHTO[node_type] * phi * fc),
        A_cs_req=C / (phi * fcu_a),
        node_type=node_type, phiFnn_aashto=phiFnn_a, node_ok_aashto=phiFnn_a >= Pn,
        phiFnn_aci=phiFnn_i, node_ok_aci=phiFnn_i >= Pn,
        eps1=strut_eps1(es, theta_deg),
        util_strut=C / phiFns_a if phiFns_a else float("inf"),
        util_node=Pn / phiFnn_a if phiFnn_a else float("inf"))
