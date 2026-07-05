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
    """使用階段混凝土壓應力上限（AASHTO 5.9.4.2）= −0.60 f'c。"""
    return -0.60 * fc


def comp_service_TW_permanent(fc: float) -> float:
    """台灣 永久載重壓應力上限 = −0.45 f'c。"""
    return -0.45 * fc


def tension_full_prestress() -> float:
    """台灣 完全預壓設計：底緣零拉 = 0.0。"""
    return 0.0


def tension_serviceIII(fc: float) -> float:
    """AASHTO Service III 底緣容許拉應力 ≈ 0.5√f'c（MPa，一般腐蝕）。"""
    return 0.5 * fc ** 0.5


# ────────────────────────────── 施拉階段（transfer）──────────────────────────────

def transfer_tension(fci: float) -> float:
    """施拉階段容許拉應力 ≈ 0.25√f'ci（MPa；AASHTO 5.9.2.3.1b 有黏結近似）。"""
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
