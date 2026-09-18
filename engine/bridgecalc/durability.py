"""耐久性設計（N1）：灌漿品質驗收、100 年設計壽命應力限制、設計年限。

對應公式卡_耐久性設計。耐久性為預防性工程，三道防線：①保護層 ②灌漿 ③應力限制控裂。
本模組將可程式化的「驗收門檻／限值」codify（值與條號寫死，比照 allowables.py）。
- 保護層數值屬查表，canonical 在「混凝土保護層快查卡」，本模組不複製。
- 混凝土疲勞壓應力（表-6.3.5）canonical 在 C1（服務性應力限制）完整表，本模組不複製（避免漂移）。
  （2026-07-01 NLM 查證道示原文：箱形 60→18.0/80→26.0；原 N1 §五 20.0/27.0 為誤，已校正與 C1 一致。）
單位：應力 MPa(=N/mm²)、強度 MPa。
"""
from dataclasses import dataclass
from typing import List

# ── 設計使用年限（年）──
DESIGN_LIFE = {"台灣": (50, 100), "AASHTO": 75, "日本": 100}

# ── 灌漿品質驗收門檻（道示Ⅲ 17.6.6(2)，三規範中最完整明文）──
GROUT = {
    "w_c_max": 0.45,            # 水灰比上限
    "f28_min": 30.0,           # 材齡 28 日壓縮強度下限 N/mm²
    "bleed_max_pct": 0.0,      # 泌水率（24h）上限 %（無泌水型）
    "vol_change_abs_pct": 0.5, # 體積變化率 ±% 以內
    "chloride_max_pct": 0.08,  # 氯離子含量（佔水泥質量）上限 %
}

# ── 100 年設計壽命鋼筋拉應力上限（道示Ⅲ 6.2.2/6.3.2）N/mm² ──
REBAR_LIMIT = {"常時": 100.0, "疲勞_一般": 180.0, "疲勞_床版翼緣": 120.0}


def design_life(code: str):
    """設計使用年限（年）。code ∈ {台灣, AASHTO, 日本}；台灣回傳 (50, 100)。"""
    return DESIGN_LIFE[code]


@dataclass
class GroutQCResult:
    all_ok: bool
    failed: List[str]   # 未達標項目名稱（空 = 全部通過）


def grout_qc_check(w_c: float, f28: float, bleed_pct: float,
                   vol_change_pct: float, chloride_pct: float) -> GroutQCResult:
    """灌漿品質驗收（道示Ⅲ 17.6.6(2)）。回傳整體判定與未達標項目清單。"""
    failed: List[str] = []
    if w_c > GROUT["w_c_max"]:
        failed.append("w/c")
    if f28 < GROUT["f28_min"]:
        failed.append("28日強度")
    if bleed_pct > GROUT["bleed_max_pct"]:
        failed.append("泌水率")
    if abs(vol_change_pct) > GROUT["vol_change_abs_pct"]:
        failed.append("體積變化")
    if chloride_pct > GROUT["chloride_max_pct"]:
        failed.append("氯離子")
    return GroutQCResult(all_ok=not failed, failed=failed)


def rebar_stress_limit(condition: str) -> float:
    """100 年耐久性鋼筋拉應力上限（道示Ⅲ 6.2.2/6.3.2）N/mm²。
    condition ∈ {常時, 疲勞_一般, 疲勞_床版翼緣}。"""
    return REBAR_LIMIT[condition]


def rebar_stress_ok(sigma_s: float, condition: str) -> bool:
    """鋼筋拉應力 σ_s 是否滿足該條件下的耐久性限值。"""
    return sigma_s <= rebar_stress_limit(condition)


def pc_fatigue_limit(Pu: float, Py: float) -> float:
    """PC 鋼材疲勞應力上限 = min(0.60·Pu, 0.75·Py)（道示Ⅲ 6.3.2 表-6.3.4）。
    Pu 極限強度、Py 降伏強度（同單位，回傳同單位）。"""
    return min(0.60 * Pu, 0.75 * Py)


# ── 台灣第十二章耐久性保護層（NLM 7d947294 核對原文 p.341–343）────────────
# 表 12.1：等級 I＝非乾濕交替（箱梁內部）；等級 II＝乾濕交替（柱、橋台、版、I/T 梁、箱梁外露面）。
# 表 12.2（一般環境）：{年限: {等級: [(最大水膠比, 最小保護層 mm), ...]}}
TW_COVER_GENERAL = {
    50: {"I": [(0.50, 30), (0.45, 25)], "II": [(0.50, 40), (0.45, 35), (0.40, 30)]},
    100: {"I": [(0.50, 35), (0.45, 30)], "II": [(0.45, 45), (0.40, 40), (0.35, 35)]},
}
# 表 12.4（鹽害）：最大水膠比、最低 f'c（kgf/cm²）
TW_SALT_MAX_WC = {"極嚴重": 0.40, "嚴重": 0.40, "中度": 0.45}
TW_SALT_MIN_FC = {"極嚴重": 350, "嚴重": 350, "中度": 280}
# 表 12.5（鹽害）主要構件最小保護層 mm：{部位: {年限: (極嚴重, 嚴重, 中度)}}
TW_COVER_SALT = {
    "基礎基樁": {50: (100, 100, 100), 100: (100, 100, 100)},
    "柱牆": {50: (100, 75, 75), 100: (100, 100, 75)},
    "橋面版頂層筋": {50: (65, 55, 50), 100: (75, 65, 60)},
    "橋面版下層筋": {50: (65, 55, 50), 100: (75, 65, 60)},
    "箱梁底層筋": {50: (65, 55, 50), 100: (75, 65, 60)},
    "梁腹版外露面": {50: (65, 55, 50), 100: (75, 65, 60)},
    "未曝露面": {50: (40, 40, 40), 100: (40, 40, 40)},
}
_SALT_IDX = {"極嚴重": 0, "嚴重": 1, "中度": 2}
PC_BASIC_COVER = 40.0   # §8.25.1：預力鋼材與主鋼筋 4 cm（基本底線）


@dataclass
class CoverReq:
    table: float        # 第十二章表列值 mm（None＝水膠比超出允許）
    basic: float        # §8.25.1 基本底線 mm
    required: float     # max(表列, 基本)
    wc_ok: bool
    source: str


def durability_cover(env: str = "general", life: int = 50, grade: str = "II",
                     wc: float = 0.45, member: str = "梁腹版外露面",
                     zone: str = "中度", basic: float = PC_BASIC_COVER) -> CoverReq:
    """台灣第十二章耐久性最小保護層（§12.4.2 表 12.2／§12.4.3 表 12.4–12.5）。

    §12.2：耐久性要求為「結構達到設計年限所需之最低要求」→ 取 max(第十二章表列, §8.25.1 基本值)。
    一般環境：表列為「最大水膠比 → 最小保護層」，實際 wc 須 ≤ 該列上限；取符合條件中最小者。
    鹽害環境：wc 須 ≤ 表 12.4 上限，保護層依表 12.5 部位×鹽害區×年限。
    ⚠ 「橋面板頂面 50 mm」是 §7.1.5 表 7.2（RC）的值，或鹽害中度 50 年的值——不是 PC 的通則。
    """
    if env == "general":
        rows = TW_COVER_GENERAL[life][grade]
        ok = [c for (w, c) in rows if wc <= w + 1e-9]
        tbl = min(ok) if ok else None
        src = f"表12.2 一般環境 等級{grade} {life}年 w/c≤{wc}"
    else:
        ok_wc = wc <= TW_SALT_MAX_WC[zone] + 1e-9
        tbl = TW_COVER_SALT[member][life][_SALT_IDX[zone]] if ok_wc else None
        ok = [tbl] if ok_wc else []
        src = f"表12.5 鹽害 {member} {zone}鹽害區 {life}年"
    req = max(tbl, basic) if tbl is not None else None
    return CoverReq(table=tbl, basic=basic, required=req, wc_ok=bool(ok), source=src)
