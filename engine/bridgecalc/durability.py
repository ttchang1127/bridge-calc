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
