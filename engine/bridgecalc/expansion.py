"""伸縮縫設計（E2）：縮短量鏈與縫口最大開度。

對應公式卡_伸縮縫設計、算例_伸縮縫設計。
最大開度 = 安裝開度 + 縮短量（溫度+潛變+乾縮；後安裝預力縮短不納入縫口計算）。
單位：mm。
"""
from dataclasses import dataclass


@dataclass
class JointResult:
    shortening: float   # 有效縮短量合計 mm
    g_max: float        # 最大開度 mm
    capacity: float     # 設計容量 mm（含裕度）
    joint_type: str     # 建議縫型


def expansion_joint(dl_T: float, dl_c: float, dl_s: float, g_install: float,
                    margin: float = 1.05) -> JointResult:
    """dl_T 溫度、dl_c 潛變、dl_s 乾縮縮短（mm）；g_install 安裝開度。"""
    shortening = dl_T + dl_c + dl_s
    g_max = g_install + shortening
    cap = round(g_max * margin)
    jt = "Strip Seal 75mm" if cap <= 75 else ("Modular ≤100mm" if cap <= 100 else "Modular >100mm")
    return JointResult(shortening, g_max, cap, jt)
