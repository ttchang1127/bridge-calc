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
    """角隅剝裂力與配筋（AASHTO Art. 5.10.9.6.5，約 0.02·P_s）。"""
    F = ratio * Ps_kN
    return SpallResult(F_spall=F, As_spall=F * 1e3 / (phi * fy))


# ── (d) 介面剪力摩擦配筋（Interface Shear Friction）────────────
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

    K1/K2 隨規範版本不同（此處預設 0.25 與 10.3 MPa 為常用值），視採用之規範版次調整。
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
