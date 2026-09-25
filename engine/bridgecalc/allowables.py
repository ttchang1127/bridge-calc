"""容許值 SSOT（值寫死＋規範出處，符合 網頁計算器SOP §八）。應力 MPa，壓為負。

本檔為全引擎「應力容許值／限值」的**單一真理源**：各設計模組的限值一律 import
本檔，不再各自硬編碼。改版只改此檔；每值標規範條號，便於送審回溯。

分類：
  使用階段 service：comp_service / comp_service_TW_permanent / tension_full_prestress / tension_serviceIII
  施拉階段 transfer：transfer_tension（拉）／transfer_comp（壓，依橋型 一般0.55/節塊0.60）
  剪力 shear：principal_tension_TW（台灣主拉）
  疲勞 fatigue：concrete_fatigue（混凝土壓疲勞）
  接縫 segmental：JOINT_MIN_COMPRESSION_MPa
  （PC 鋼材／鋼筋疲勞、灌漿驗收等 codified 限值見 durability.py 之 REBAR_LIMIT/GROUT/DESIGN_LIFE，
   已按條號整理，本檔以引用方式對齊，不重複定義。）
"""
import math

# ────────────────────────────── 使用階段（service）──────────────────────────────

def comp_service(fc: float) -> float:
    """使用階段混凝土壓應力上限（AASHTO 5.9.4.2；台灣 §8.15.2 2.(1) 同值）= −0.60 f'c。

    台灣原文：「(1) 除本項(2)、(3)所述之情形外，所有載重組合下壓應力應小於 0.6f'c。」
    """
    return -0.60 * fc


def comp_service_TW_permanent(fc: float) -> float:
    """台灣 有效預力＋永久靜載重之壓應力上限〔§8.15.2 2.(2)，p.155–156〕= −0.40 f'c。

    🔴 **2026-09-24 稽核更正**：原值寫 −0.45 f'c 且未標條號。0.45 f'c 是 **AASHTO**
      5.9.4.2.1 的「有效預力＋永久載重」值；台灣 §8.15.2 2.(2) 原文為
      「有效預力加上永久靜載重產生之壓應力應小於 **0.40f'c**」——兩制不同值。
      誤用 AASHTO 值會放寬 12.5%（f'c=40 時 18.0 vs 16.0 MPa），**偏不安全**。
      （本函式在更正時尚無呼叫者，故 golden 不受影響；決策 19 同型誤植。）
    """
    return -0.40 * fc


def comp_service_TW_live_half(fc: float) -> float:
    """台灣 §8.15.2 2.(3)：「(2)所得計算壓應力之半加上活重產生之壓應力應小於 0.40f'c」。

    ⚠ 這是台灣**獨有的第三道**檢核（AASHTO 無對應項），先前引擎完全沒有。
    檢核式：0.5·σ_(永久) + σ_(活載) ≤ 0.40 f'c（皆取壓應力絕對值）。
    """
    return -0.40 * fc


def tension_precompressed_TW(fc: float, bonded: bool = True,
                             segmental: bool = False, corrosive: bool = False) -> float:
    """台灣 §8.15.2 2. 預壓拉力區之容許拉應力（損失後，MPa，拉為正）。

    原文（kgf/cm² 制，括號為 MPa 制）：
      配置握裹鋼筋(一般橋梁) 1.6√f'c (0.498√f'c)；配置握裹鋼筋(節塊橋梁) 0.8√f'c (0.249√f'c)；
      曝露於嚴重腐蝕情況，如沿海地區者(一般橋梁) 0.8√f'c (0.249√f'c)；(節塊橋梁) 0；
      不配置握裹鋼筋之構材 0。
    💡 一般橋梁有握裹之 0.498√f'c 與 AASHTO Service III 的 0.5√f'c **幾乎同值**（決策 23 之跨軌吻合）。
    """
    if not bonded:
        return 0.0
    if corrosive:
        return 0.0 if segmental else 0.249 * math.sqrt(fc)
    return 0.249 * math.sqrt(fc) if segmental else 0.498 * math.sqrt(fc)


# 台灣 §8.15.2 1. 拉力之絕對上限（「14 kgf/cm²」）換算 MPa
TW_TRANSFER_TENSION_CAP_MPa = 14.0 * 0.0980665     # = 1.3729 MPa


def transfer_tension_TW(fci: float, bonded: bool = False,
                        segmental: bool = False) -> float:
    """台灣 §8.15.2 1. 拉力（損失前暫時應力，「其他區域」），MPa（拉為正）。

    原文：不配置握裹鋼筋之拉力區　一般橋梁「**14 kgf/cm² 或 0.8√f'ci**」／節塊橋梁「不容許產拉應力」；
          配置握裹鋼筋之拉力區　一般橋梁 2.0√f'ci (0.623√f'ci)／節塊橋梁 1.6√f'ci (0.498√f'ci)。
    🔴 無握裹者為 **min(0.25√f'ci, 1.373 MPa)**——0.8√f'ci〔kgf〕換算即 0.2505√f'ci〔MPa〕，
      與 AASHTO 5.9.4.1.2 的 min(0.25√f'ci, 1.38 MPa) 實質同值；**兩制都有上限，`transfer_tension`
      卻沒有**（f'ci=32 時 1.414 vs 1.373，寬鬆 3.0%）。詳見該函式之註記。
    """
    if segmental:
        return 0.0 if not bonded else 0.498 * math.sqrt(fci)
    if bonded:
        return 0.623 * math.sqrt(fci)
    return min(0.25 * math.sqrt(fci), TW_TRANSFER_TENSION_CAP_MPa)


def tension_full_prestress() -> float:
    """台灣 完全預壓設計：底緣零拉 = 0.0。"""
    return 0.0


def tension_serviceIII(fc: float) -> float:
    """AASHTO Service III 底緣容許拉應力 ≈ 0.5√f'c（MPa，一般腐蝕）。"""
    return 0.5 * fc ** 0.5


# ────────────────────────────── 施拉階段（transfer）──────────────────────────────

def transfer_tension(fci: float) -> float:
    """施拉階段容許拉應力 ≈ 0.25√f'ci（MPa；AASHTO 5.9.2.3.1b 有黏結近似）。

    ⚠ **2026-09-24 稽核發現：本式缺上限。** AASHTO 5.9.4.1.2 與台灣 §8.15.2 1. 對「其他區域·
      無握裹鋼筋」都是 **min(0.25√f'ci, ≈1.38 MPa)**（台灣原文「14 kgf/cm² 或 0.8√f'ci」，
      14 kgf/cm²＝1.3729 MPa）。f'ci > 30.1 MPa 時上限才控制——40m 參考橋 f'ci=32 恰落在
      上限側（1.414 vs 1.373，本式寬鬆 3.0%、偏不安全）。
      本檔 `neg_top_tension_limit`（staging.py）同一限值**有**帶上限，引擎內部不一致。
      🔴 加上限會動既有 golden（施拉階段判定），故保留原式待裁示；
      需要台灣明文值者請用 `transfer_tension_TW()`。
    """
    return 0.25 * math.sqrt(fci)


# 施拉階段容許壓應力係數（公式卡_服務性應力限制 C1「施拉壓應力」表；後拉法）
TRANSFER_COMP_FACTOR = {"一般": 0.55, "節塊": 0.60}


def transfer_comp(fci: float, bridge_type: str = "一般") -> float:
    """施拉階段容許壓應力（C1，壓為負）。一般橋梁 −0.55 f'ci／節塊橋梁 −0.60 f'ci。

    ⚠️ 40m 參考橋為場鑄一般橋梁，C1 定為 0.55 f'ci；惟現行施工模組（construction.py）
    與其算例／golden 沿用 0.60（節塊值），S2 底板「−19.18 剛好過 fci32(19.2)」即依 0.60。
    此為待決工程判斷（見待辦），故本函式以 bridge_type 參數化，預設 0.55。
    """
    return -TRANSFER_COMP_FACTOR[bridge_type] * fci


# ────────────────────────────── 撓度（deflection）──────────────────────────────
# 台灣同一組限值在三章各有條文、**數值相同但條號不同**（2026-09-24 NLM 核）：
#   RC §7.3.12（p.141）／PC §8.11.3（p.152–153）／鋼結構 §9.1.7（p.194，語氣為「以…為宜」）
# 分母 canonical：簡支/連續 800；市區部分供人行 1000；懸臂端 300；懸臂供人行道 375。
TW_DEFLECTION_DENOM = {"一般": 800, "市區人行": 1000, "懸臂": 300, "懸臂人行": 375}


def deflection_limit_TW(L: float, case: str = "一般") -> float:
    """台灣活載重＋衝擊力之撓度限值（mm，L 同單位）。

    case：'一般'（簡支或連續，L/800）／'市區人行'（L/1000）／
          '懸臂'（臂長 L_c/300）／'懸臂人行'（L_c/375）。
    ⚠ 懸臂者 L 請傳**懸臂長度 L_c**，不是跨徑。
    """
    if case not in TW_DEFLECTION_DENOM:
        raise ValueError(f"case 須為 {list(TW_DEFLECTION_DENOM)}")
    return L / TW_DEFLECTION_DENOM[case]


def long_term_factor_TW(by_Ig: bool = True, As_comp_ratio: float = 0.0) -> float:
    """台灣 RC 長期載重因素〔§7.1.22 7.(4)，p.111〕。

    原文：「a. 如即時撓度係按全斷面慣性矩 Ig 求得時，長期載重因素可取 4。
            b. 如即時撓度係按有效慣性矩 Ie 求得時，長期載重因素可取 3−1.2(A's/As) ≥ 1.6」
    🔴 **下限 1.6 不可漏**：A's/As > 1.17 時 3−1.2r 即低於 1.6，漏掉下限會低估長期撓度（偏不安全）。
    ⚠ 本式屬**第七章 RC**；預力構材之長期撓度請用潛變係數逐步計算（§8.11.1 明列應考慮之影響項）。
    """
    return 4.0 if by_Ig else max(3.0 - 1.2 * As_comp_ratio, 1.6)


# ────────────────────────────── 剪力（shear）──────────────────────────────

def principal_tension_TW(fc: float) -> float:
    """台灣主拉應力容許 = 0.3√f'c(kgf/cm²)，換算 MPa = 0.094√f'c。"""
    return 0.094 * math.sqrt(fc)


# ────────────────────────────── 疲勞（fatigue）──────────────────────────────

def concrete_fatigue(fc: float) -> float:
    """混凝土疲勞壓應力上限（AASHTO 5.5.3.1）= 0.40 f'c（正值，與 σ_c,max 絕對值比）。"""
    return 0.40 * fc


# ────────────────────────────── 接縫（segmental）──────────────────────────────

# 節塊接縫拼裝期最小壓應力（AASHTO/道示 ≈ 30 psi）
JOINT_MIN_COMPRESSION_MPa = 0.21
