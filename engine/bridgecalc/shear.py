"""腹板抗剪（D1）：軸壓 fpc、鋼腱垂直分量 Vp、莫耳圓主拉應力、Vcw、所需箍筋。

對應公式卡 D1（腹板剪力設計）、算例_腹板抗剪設計。單位：力 N、長度 mm、應力 MPa。
台灣主拉/Vcw 係數 0.3√f'c(kgf/cm²) 已換算為 MPa：0.094·√f'c（f'c MPa）。
"""
from dataclasses import dataclass
from math import sqrt
from .model import Section
from . import allowables


@dataclass
class ShearResult:
    fpc: float        # 軸向壓應力 MPa（正值）
    slope: float      # 控制斷面鋼腱斜率（≈sinθ）
    Vp: float         # 鋼腱垂直分量（每腹板）N
    tau: float        # 剪應力 MPa
    sigma1: float     # 主拉應力 MPa
    sigma1_limit: float   # 台灣主拉容許 MPa
    sigma1_ok: bool       # σ1 ≤ 限值？（近支承常超→靠箍筋）
    Vcw: float        # 混凝土抗剪（含 Vp）N
    Vs_req: float     # 所需鋼材抗剪 N
    Av_s_req: float   # 所需箍筋 mm²/mm


def principal_tension_limit_TW(fc: float) -> float:
    """台灣主拉應力容許（→ allowables.principal_tension_TW，0.094√f'c）。"""
    return allowables.principal_tension_TW(fc)


def shear_web(Pe: float, section: Section, e: float, fc: float,
              bw_eff: float, dv: float, Vu: float,
              x_control: float, L: float,
              n_webs: int = 2, phi: float = 0.85, fsy: float = 420.0) -> ShearResult:
    """近支承控制斷面腹板抗剪（每腹板）。

    Pe 全斷面 N、e 跨中偏心 mm、Vu 每腹板設計剪力 N、x_control 控制斷面距支承 mm。
    """
    slope = (4 * e / L) * (1 - 2 * x_control / L)        # 拋物線斜率 ≈ sinθ
    return shear_web_at(Pe, section, slope, fc, bw_eff, dv, Vu, n_webs, phi, fsy)


def shear_web_at(Pe: float, section: Section, slope: float, fc: float,
                 bw_eff: float, dv: float, Vu: float,
                 n_webs: int = 2, phi: float = 0.85, fsy: float = 420.0) -> ShearResult:
    """任意斷面腹板抗剪（每腹板；連續梁墩兩側等）。

    slope＝鋼腱斜率 de/dx（e 形心下為正、x 沿橋，無因次）；Vu＝每腹板設計剪力 N（帶號，dM/dx 慣例）。
    預力垂直分力 Vp = (Pe/n)·|slope|，**僅當與 Vu 反向（即 slope 與 Vu 同號）時計為有利**，否則取負值（不利）。
    簡支近左支承：e 向跨中增加 → slope>0、Vu>0 → 有利（同 shear_web）。
    連續梁墩左側：腱往墩頂上升 → slope<0、Vu<0 → 有利。
    """
    fpc = Pe / section.A
    Vabs = abs(Vu)
    Vp = (Pe / n_webs) * abs(slope) * (1 if slope * Vu >= 0 else -1)
    tau = Vabs / (bw_eff * dv)
    sigma1 = -fpc / 2 + sqrt((fpc / 2) ** 2 + tau ** 2)
    lim = principal_tension_limit_TW(fc)
    vc = 0.094 * sqrt(fc)                                 # 台灣混凝土拉貢獻 MPa
    Vcw = (vc + 0.3 * fpc) * bw_eff * dv + Vp
    Vs_req = (Vabs - phi * Vcw) / phi
    Av_s_req = max(Vs_req, 0.0) / (fsy * dv)
    return ShearResult(fpc, slope, Vp, tau, sigma1, lim, sigma1 <= lim,
                       Vcw, Vs_req, Av_s_req)


def phiVn(Vcw: float, Av_s_provided: float, dv: float,
          phi: float = 0.85, fsy: float = 420.0) -> float:
    """實配箍筋下的設計抗剪 φ(Vc+Vs)。"""
    Vs = Av_s_provided * fsy * dv
    return phi * (Vcw + Vs)


def Av_s_min_TW(fc: float = None, bw_eff: float = 0.0, fsy: float = 420.0) -> float:
    """台灣最小腹板（剪力）鋼筋 Av/s = 0.345·b′/fsy（mm²/mm）。

    原文（2026-09-20 NLM 核）：§8.20.3 3.「腹板鋼筋最小斷面積應為 Av = 3.5·b′s/fsy
    （b′、s 單位 cm，fsy 單位 kgf/cm²）＝ 0.345·b′s/fsy（MPa 制）」（第八章 p.167，式 8-31）；
    第七章 §7.1.9 1.(2) RC 同式。**兩章均無含 √f'c 之項**。
    fc 保留為相容參數（不影響結果），台灣式與混凝土強度無關。

    🔴 2026-09-20 更正：原式 max(0.2√f'c·b/fsy, 0.35·b/fsy) 有兩個錯——①√f'c 項非台灣規範
    （出自 ACI 318／AASHTO）；②0.2 是 **kgf/cm² 制**係數卻餵 MPa（ACI SI 制為 0.062、
    AASHTO 5.8.2.5 為 0.083），f'c=40 時要求量為明文值的 3.7 倍。
    使用者裁示（2026-09-20）：判定用台灣明文，AASHTO 以 Av_s_min_AASHTO 並列參考。
    """
    return 0.345 * bw_eff / fsy


def Av_s_min_AASHTO(fc: float, bw_eff: float, fy: float = 420.0) -> float:
    """AASHTO LRFD 5.8.2.5 最小剪力鋼筋 Av/s = 0.083√f'c·b_v/f_y（mm²/mm）。參考用，非台灣規範。"""
    return 0.083 * sqrt(fc) * bw_eff / fy


STD_STIRRUP_SPACINGS = (100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 600)


def stirrup_max_spacing_TW(Vs: float, fc: float, bw: float, dv: float, h: float):
    """台灣橋規 §8.20.3 箍筋最大間距與斷面上限（檢核流程_腹板抗剪 STEP 6）。

    Vs 每腹板所需 N。回傳 (s_max mm, 是否減半, 斷面足夠 Vs ≤ 0.66√f'c·b·d_v)。
    V_s ≤ 0.33√f'c·b·d_v → min(0.75h, 600)；超過 → min(0.375h, 300)。
    """
    lim1 = 0.33 * sqrt(fc) * bw * dv
    lim2 = 0.66 * sqrt(fc) * bw * dv
    halved = Vs > lim1
    s_max = min(0.375 * h, 300.0) if halved else min(0.75 * h, 600.0)
    return s_max, halved, Vs <= lim2


def stirrup_pick_spacing(s_allow: float, spacings=STD_STIRRUP_SPACINGS):
    """取不大於 s_allow 的最大標準間距；皆不滿足回傳 None（須改箍筋號數／肢數）。"""
    ok = [s for s in spacings if s <= s_allow + 1e-9]
    return max(ok) if ok else None


def stirrup_zones(xs: list, picks: list, supports: list = None):
    """由逐斷面選定間距產生分區：相鄰取樣點之間取兩端較密者，連續同間距合併。

    xs 由小到大；picks 對應間距（None＝不足）；supports 為支承里程——支承面至第一個取樣點
    沿用該取樣點間距（距支承 d_v 內以 d_v 斷面設計）。回傳 [(x1, x2, s), ...]。
    """
    segs = []
    sup = list(supports or [])
    for i in range(len(xs) - 1):
        if any(xs[i] < sx < xs[i + 1] for sx in sup):
            continue                                   # 跨過支承的區間改由下方支承面處理
        a, b = picks[i], picks[i + 1]
        s = None if (a is None or b is None) else min(a, b)
        segs.append([xs[i], xs[i + 1], s])
    if supports:
        for sx in supports:
            left = [i for i, x in enumerate(xs) if x < sx]
            right = [i for i, x in enumerate(xs) if x > sx]
            if right and (not left or xs[left[-1]] < sx):
                j = right[0]
                segs.append([sx, xs[j], picks[j]])
            if left and (not right or xs[right[0]] > sx):
                j = left[-1]
                segs.append([xs[j], sx, picks[j]])
    segs.sort(key=lambda t: (t[0], t[1]))
    merged = []
    for a, b, s in segs:
        if b - a < 1e-9:
            continue
        if merged and abs(merged[-1][1] - a) < 1e-9 and merged[-1][2] == s:
            merged[-1][1] = b
        else:
            merged.append([a, b, s])
    return [tuple(m) for m in merged]
