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


def _max_moving(L: float, il, axles, spacing, step: float) -> float:
    """軸組沿梁掃描取 |ΣPᵢ·η| 最大者（保留正負號）。

    車可雙向行駛 → 原向與掉頭（軸序與軸距鏡射）都掃；單向會使非跨中斷面
    彎矩與短跨支承剪力偏低。整數步進免浮點累加漂移；軸位夾到 [0, L]，
    否則壓在支承上的軸會因 −1e-15 被略去。
    """
    span = spacing[-1]
    mirrored = (tuple(reversed(axles)), tuple(span - d for d in reversed(spacing)))
    n = int(round((L + span) / step))
    best = 0.0
    for P_set, d_set in ((tuple(axles), tuple(spacing)), mirrored):
        for k in range(n + 1):
            s = -span + k * step
            tot = sum(P * il(min(max(s + d, 0.0), L))
                      for P, d in zip(P_set, d_set) if -1e-9 <= s + d <= L + 1e-9)
            if abs(tot) > abs(best):
                best = tot
    return best


def max_moment_moving(L: float, a: float, axles=HL93_AXLES,
                      spacing=HL93_SPACING, step: float = 0.1) -> float:
    """移動車組於斷面 a 之最大彎矩 max ΣPᵢ·η（雙向掃描）。"""
    return _max_moving(L, lambda p: il_moment_simple(L, a, p), axles, spacing, step)


def abs_max_moment(L: float, axles=HL93_AXLES, spacing=HL93_SPACING,
                   step: float = 0.2) -> float:
    """絕對最大彎矩：掃描斷面 a 與車組位置，取全梁最大。"""
    best = 0.0
    j = 1
    while j * step < L - 1e-9:
        m = max_moment_moving(L, j * step, axles, spacing, step)
        if abs(m) > abs(best):
            best = m
        j += 1
    return best


def lane_moment_simple(L: float, w: float = HL93_LANE) -> float:
    """車道均布載重跨中彎矩 = w·L²/8（= w × 影響線正面積）。"""
    return w * L**2 / 8


def moment_envelope_simple(L: float, axles=HL93_AXLES, spacing=HL93_SPACING,
                           n: int = 20):
    """簡支梁設計卡車彎矩包絡線：每斷面 a 的最大彎矩。

    回傳 [(a, M_max), ...]（n+1 點）。峰值 = 絕對最大彎矩（≈跨中）。
    """
    return [(L * i / n, max_moment_moving(L, L * i / n, axles, spacing))
            for i in range(n + 1)]


def hl93_per_lane_moment(L: float) -> float:
    """每設計車道 M_LL+IM = (1+IM)·卡車絕對最大 + 車道（簡支跨中近似）。"""
    return (1 + IM) * abs_max_moment(L) + lane_moment_simple(L)


# ── 台灣 HS20-44（公路橋梁設計規範 §3.6/§3.13）──
# 設計車輛 + 車道載重「取大」（非 HL-93 並計）；衝擊及於控制者。
TW_HS20_AXLES = (36.0, 144.0, 144.0)
TW_HS20_SPACING = (0.0, 4.25, 8.5)
TW_LANE_W = 9.4       # 車道均布載重 kN/m
TW_LANE_PM = 80.0     # 車道集中載重（彎矩用）kN
TW_LANE_PV = 116.0    # 車道集中載重（剪力用）kN


def taiwan_impact(L: float) -> float:
    """台灣衝擊係數 I = 15.24/(L+38.1) ≤ 0.30（§3.13）。"""
    return min(0.30, 15.24 / (L + 38.1))


def max_shear_moving(L: float, a: float, axles=TW_HS20_AXLES,
                     spacing=TW_HS20_SPACING, step: float = 0.1) -> float:
    """移動車組於斷面 a 之最大剪力 max ΣPᵢ·η_V（雙向掃描）。a=0 即支承反力。"""
    return _max_moving(L, lambda p: il_shear_simple(L, a, p), axles, spacing, step)


def taiwan_lane_moment(L: float) -> float:
    """台灣車道載重跨中彎矩 = w·L²/8 + 集中載重 P_M·L/4。"""
    return TW_LANE_W * L**2 / 8 + TW_LANE_PM * L / 4


def taiwan_lane_shear(L: float) -> float:
    """台灣車道載重支承剪力 = w·L/2 + 集中載重 P_V。"""
    return TW_LANE_W * L / 2 + TW_LANE_PV


def taiwan_truck_moment(L: float) -> float:
    """台灣 HS20-44 設計卡車絕對最大彎矩（不含衝擊）。"""
    return abs_max_moment(L, TW_HS20_AXLES, TW_HS20_SPACING)


def taiwan_moment_envelope(L: float, n: int = 20):
    """HS20-44 設計卡車彎矩包絡線 [(a, M_max), ...]。"""
    return moment_envelope_simple(L, TW_HS20_AXLES, TW_HS20_SPACING, n)


def taiwan_per_lane_moment(L: float) -> float:
    """台灣每設計車道 M_LL+IM = max(設計卡車絕對最大, 車道載重) × (1+I)。"""
    truck = abs_max_moment(L, TW_HS20_AXLES, TW_HS20_SPACING)
    return max(truck, taiwan_lane_moment(L)) * (1 + taiwan_impact(L))


def taiwan_truck_shear(L: float) -> float:
    """台灣 HS20-44 設計卡車支承最大剪力（不含衝擊）。L ≥ 8.5 m 時 = 324 − 918/L（後軸壓支承、前軸朝跨內）。"""
    return abs(max_shear_moving(L, 0.0, TW_HS20_AXLES, TW_HS20_SPACING))


def taiwan_per_lane_shear(L: float) -> float:
    """台灣每設計車道 V_LL+IM = max(設計卡車支承剪力, 車道剪力) × (1+I)。"""
    return max(taiwan_truck_shear(L), taiwan_lane_shear(L)) * (1 + taiwan_impact(L))
