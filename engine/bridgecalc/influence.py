"""影響線（簡支梁解析）與移動車組最大效應。

對應快查卡_影響線方法、算例_40m參考橋活載影響線實算。
JS 網頁計算器（數值剛度法）以此模組的解析值對齊（簡支三角形/階梯）。
單位：長度 m、載重 kN → 彎矩 kN·m、剪力 kN。
"""
from dataclasses import dataclass

# HL-93 設計卡車（軸載 kN、相對軸距 m；求最大彎矩取最小軸距 4.3）
HL93_AXLES = (35.0, 145.0, 145.0)
HL93_SPACING = (0.0, 4.3, 8.6)
HL93_LANE = 9.3      # 車道載重 kN/m
IM = 0.33            # 衝擊（卡車，非車道）


def il_moment_simple(L: float, a: float, p: float) -> float:
    """簡支梁斷面 a 之彎矩影響線縱距（單位載重於 p）。三角形，峰值 a(L−a)/L。"""
    if p <= a:
        return p * (L - a) / L
    return a * (L - p) / L


def il_shear_simple(L: float, a: float, p: float) -> float:
    """剪力影響線縱距（a 右側為正）。a 處跳躍 1。"""
    if p < a:
        return -p / L
    return (L - p) / L


def il_moment_peak(L: float, a: float) -> float:
    """彎矩影響線峰值 = a(L−a)/L（載重在 a）。"""
    return a * (L - a) / L


def max_moment_moving(L: float, a: float, axles=HL93_AXLES,
                      spacing=HL93_SPACING, step: float = 0.1) -> float:
    """移動車組於斷面 a 之最大彎矩 max ΣPᵢ·η（卡車沿梁掃描）。"""
    best = 0.0
    s = -spacing[-1]
    while s <= L:
        tot = sum(P * il_moment_simple(L, a, s + dx)
                  for P, dx in zip(axles, spacing) if 0 <= s + dx <= L)
        if abs(tot) > abs(best):
            best = tot
        s += step
    return best


def abs_max_moment(L: float, axles=HL93_AXLES, spacing=HL93_SPACING,
                   step: float = 0.2) -> float:
    """絕對最大彎矩：掃描斷面 a 與車組位置，取全梁最大。"""
    best = 0.0
    a = step
    while a < L:
        m = max_moment_moving(L, a, axles, spacing, step)
        if abs(m) > abs(best):
            best = m
        a += step
    return best


def lane_moment_simple(L: float, w: float = HL93_LANE) -> float:
    """車道均布載重跨中彎矩 = w·L²/8（= w × 影響線正面積）。"""
    return w * L**2 / 8


def hl93_per_lane_moment(L: float) -> float:
    """每設計車道 M_LL+IM = (1+IM)·卡車絕對最大 + 車道（簡支跨中近似）。"""
    return (1 + IM) * abs_max_moment(L) + lane_moment_simple(L)
