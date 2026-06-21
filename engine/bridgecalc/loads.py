"""活載與 AASHTO 載重組合（L0）。彎矩單位 kN·m。

對應：算例_40m參考橋活載基準統一、算例_40m參考橋載重組合。
"""


def lane_live_load(per_lane_M: float, n_lanes: int, m: float = 1.0) -> float:
    """活載彎矩 = 每車道 M_LL+IM × 設計車道數 × 多車道係數 m。

    per_lane_M 為每設計車道含衝擊 IM 的彎矩（影響線實算，5,673 kN·m）。
    2 設計車道 m=1.0 → 11,346 kN·m（全箱）。
    """
    return per_lane_M * n_lanes * m


def combinations(M_DC: float, M_DW: float, M_LL_IM: float) -> dict:
    """AASHTO §3.4.1 載重組合（簡支，無次彎矩 PS）。"""
    return {
        "Strength_I":  1.25 * M_DC + 1.50 * M_DW + 1.75 * M_LL_IM,
        "Service_I":   M_DC + M_DW + M_LL_IM,
        "Service_III": M_DC + M_DW + 0.8 * M_LL_IM,   # PC 拉應力
    }
