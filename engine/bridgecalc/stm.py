"""壓拉桿模型設計（F2）：錨碇區 General Zone 爆裂力、壓桿/拉桿/節點承載力。

對應公式卡_STM壓拉桿模型設計、算例_STM端橫隔版設計。
General Zone 爆裂力（AASHTO 5.9.5.6.3）：T_burst = 0.25·ΣP·(1−a/h)，
分佈深度 d_burst = 0.5(h−2e)。壓桿/節點以有效抗壓強度 f_cu = 0.85·β·f'c 驗核：
  β_strut：0.60 無束制 / 0.75 有橫向拉力配筋 / 1.00 橫向束制；
  β_node ：CCC 1.00 / CCT 0.80 / CTT 0.60 / TTT 0.40。
拉桿 A_s = T/(φ·f_y)，φ=0.90；壓桿與節點 φ=0.70。
單位：力 N、長度 mm、應力 MPa。
"""
from dataclasses import dataclass

PHI_STM = 0.70      # 壓桿/節點折減係數
PHI_TIE = 0.90      # 拉桿（鋼筋）折減係數
BETA_NODE = {"CCC": 1.00, "CCT": 0.80, "CTT": 0.60, "TTT": 0.40}
BETA_STRUT = {"unconfined": 0.60, "reinforced": 0.75, "confined": 1.00}


def burst_force(P: float, a: float, h: float) -> float:
    """General Zone 爆裂力 T_burst = 0.25·ΣP·(1 − a/h)。P[N] 錨座合力、a 錨座分佈寬、h 斷面高。"""
    return 0.25 * P * (1 - a / h)


def burst_depth(h: float, e_anc: float) -> float:
    """爆裂力分佈中心深度 d_burst = 0.5·(h − 2·e_anc)，mm。"""
    return 0.5 * (h - 2 * e_anc)


def f_cu(fc: float, beta: float) -> float:
    """有效抗壓強度 f_cu = 0.85·β·f'c，MPa。"""
    return 0.85 * beta * fc


def strut_capacity(fc: float, beta_s: float, A_cs: float, phi: float = PHI_STM) -> float:
    """壓桿承載力 φF_ns = φ·f_cu·A_cs，N。A_cs[mm²]。"""
    return phi * f_cu(fc, beta_s) * A_cs


def node_capacity(fc: float, node_type: str, A_n: float, phi: float = PHI_STM) -> float:
    """節點承載力 φF_nn = φ·0.85·β_n·f'c·A_n，N。node_type ∈ {CCC,CCT,CTT,TTT}。"""
    return phi * f_cu(fc, BETA_NODE[node_type]) * A_n


def tie_reinforcement(T: float, fy: float = 420, phi: float = PHI_TIE) -> float:
    """拉桿鋼筋需求 A_s = T / (φ·f_y)，mm²。"""
    return T / (phi * fy)


@dataclass
class STMResult:
    sigma_pe: float        # 遠端均佈預壓 Pe/Ac，MPa
    T_burst: float         # General Zone 爆裂力 N
    d_burst: float         # 爆裂力分佈深度 mm
    As_burst: float        # 爆裂鋼筋需求 mm²
    strut_force: float     # 主壓桿力 N
    beta_strut_required: str  # 主壓桿所需 β 等級
    strut_ok: bool         # 主壓桿（採所需 β）是否通過


def general_zone_burst(P_web: float, a: float, h: float, e_anc: float,
                       Pe: float, Ac: float, A_cs_strut: float,
                       strut_force: float, fc: float = 40, fy: float = 420) -> STMResult:
    """端橫隔版 General Zone STM 驗算。

    P_web[N]：單腹板錨座合力；a[mm]：錨座分佈寬；h[mm]：斷面高；e_anc[mm]：錨座群形心偏心；
    Pe[N]/Ac[mm²]：遠端均佈預壓；A_cs_strut[mm²]：主壓桿有效斷面；strut_force[N]：主壓桿力。
    主壓桿 β 自 0.60→0.75→1.00 逐級提升，取首個通過者；皆不過則回報 confined 不足。
    """
    Tb = burst_force(P_web, a, h)
    db = burst_depth(h, e_anc)
    As = tie_reinforcement(Tb, fy)
    sig = Pe / Ac
    beta_req, ok = "confined", False
    for name in ("unconfined", "reinforced", "confined"):
        if strut_capacity(fc, BETA_STRUT[name], A_cs_strut) >= strut_force:
            beta_req, ok = name, True
            break
    return STMResult(sigma_pe=sig, T_burst=Tb, d_burst=db, As_burst=As,
                     strut_force=strut_force, beta_strut_required=beta_req, strut_ok=ok)


# ── AASHTO LRFD 2008 SI 原文版（2026-09-21 NLM 核）──────────────────────────
# 🔴 本檔上方的 BETA_STRUT／BETA_NODE（0.60/0.75/1.00、1.00/0.80/0.60/0.40）是 **ACI 318 式**
#   的 β 分類表。**AASHTO LRFD 2008 SI 沒有 β_s 表**——壓桿以應變相容式軟化：
#     f_cu = f'c / (0.8 + 170·ε₁) ≤ 0.85 f'c                      （Eq. 5.6.3.3.3-1，p.5-32）
#     ε₁ = ε_s + (ε_s + 0.002)·cot²α_s                            （Eq. 5.6.3.3.3-2）
#       α_s＝壓桿與相鄰拉桿之最小夾角；ε_s＝拉桿方向之混凝土拉應變。
#   節點區則直接給容許壓應力（Art. 5.6.3.5，p.5-34）：
#     CCC 0.85φf'c／CCT 0.75φf'c／CTT 0.65φf'c    （φ＝0.70，Art. 5.5.4.2.1 p.5-27~28）
#   換算成「0.85β」形式的等效 β 為 CCC 1.00／CCT 0.88／CTT 0.76——**與 ACI 的 0.80／0.60 不同**。
#   拉桿：P_n = f_y·A_st + A_ps(f_pe + f_y)（Eq. 5.6.3.4.1-1，p.5-33；無一般鋼筋時 f_y 可取 420）。
NODE_LIMIT_AASHTO = {"CCC": 0.85, "CCT": 0.75, "CTT": 0.65}
EPS_S_YIELD = 420.0 / 200000.0      # 拉桿降伏時之鋼筋應變 0.0021（f_y/E_s）


def strut_eps1(eps_s: float, alpha_s_deg: float) -> float:
    """主拉應變 ε₁ = ε_s + (ε_s + 0.002)·cot²α_s（Eq. 5.6.3.3.3-2）。"""
    import math as _m
    a = _m.radians(alpha_s_deg)
    return eps_s + (eps_s + 0.002) / (_m.tan(a) ** 2)


def strut_fcu_aashto(fc: float, eps_s: float = EPS_S_YIELD,
                     alpha_s_deg: float = 45.0) -> float:
    """壓桿有效抗壓強度 f_cu = f'c/(0.8+170ε₁) ≤ 0.85f'c（Eq. 5.6.3.3.3-1），MPa。

    ⚠ α_s 越小（壓桿與拉桿越接近平行）→ ε₁ 越大 → f_cu 越低；α_s = 90° 時退化為 f'c/(0.8+170ε_s)。
    """
    e1 = strut_eps1(eps_s, alpha_s_deg)
    return min(fc / (0.8 + 170.0 * e1), 0.85 * fc)


def strut_capacity_aashto(fc: float, A_cs: float, eps_s: float = EPS_S_YIELD,
                          alpha_s_deg: float = 45.0, phi: float = PHI_STM) -> float:
    """壓桿承載力 φF_ns = φ·f_cu·A_cs（AASHTO 2008 式），N。"""
    return phi * strut_fcu_aashto(fc, eps_s, alpha_s_deg) * A_cs


def node_capacity_aashto(fc: float, node_type: str, A_n: float,
                         phi: float = PHI_STM) -> float:
    """節點承載力 φF_nn＝{0.85／0.75／0.65}·φ·f'c·A_n（Art. 5.6.3.5），N。node_type∈{CCC,CCT,CTT}。"""
    return NODE_LIMIT_AASHTO[node_type] * phi * fc * A_n


def tie_strength_aashto(A_st: float, fy: float = 420.0, A_ps: float = 0.0,
                        f_pe: float = 0.0) -> float:
    """拉桿公稱拉力 P_n = f_y·A_st + A_ps(f_pe + f_y)（Eq. 5.6.3.4.1-1），N。"""
    return fy * A_st + A_ps * (f_pe + fy)
