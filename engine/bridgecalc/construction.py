"""施工階段應力（H1 全支撐 / H2 MSS 移動模架）：支架上施拉的頂緣引張支配條件。

對應算例_40m參考橋施工階段應力歷程、公式卡_全支撐工法/移動模架工法、C1 §六。
關鍵：梁在支架/托架上施拉時，自重彎矩尚未活化（M_sw≈0，由支架承擔），預力大偏心
使頂緣淨受拉（過平衡 LBR>1）→ 須分批張拉降低瞬時預力；脫架後自重活化 → 斷面回壓。
複用 service.stresses()（同符號慣例，壓為負）。單位：力 N、彎矩 kN·m、應力 MPa。
"""
import math
from dataclasses import dataclass

from .model import Section
from .service import stresses


@dataclass
class StageStress:
    sigma_top: float    # 頂緣應力 MPa（壓為負）
    sigma_bot: float    # 底緣應力 MPa
    top_ok: bool        # 頂緣拉應力 ≤ 施拉容許拉
    bot_ok: bool        # 底緣壓應力 ≥ 施拉容許壓（壓不超限）


def transfer_tension_limit(fci: float) -> float:
    """施拉階段容許拉應力 ≈ 0.25·√f'ci（MPa；AASHTO 5.9.2.3.1b 有黏結近似）。"""
    return 0.25 * math.sqrt(fci)


def transfer_comp_limit(fci: float) -> float:
    """施拉階段容許壓應力 −0.60·f'ci（MPa，壓為負）。"""
    return -0.60 * fci


def stage_stress(P: float, sec: Section, e: float, M_sw_kNm: float,
                 fci: float) -> StageStress:
    """施工階段跨中頂/底緣應力與判定。

    P[N]：本階段有效預力（分批時 = Pi × 已張組數比例）；
    M_sw_kNm：該階段「已活化」的自重彎矩（支架/托架上 ≈ 0；脫架後 = M_DC）。
    """
    st, sb = stresses(P, sec, e, M_sw_kNm)
    return StageStress(
        sigma_top=st, sigma_bot=sb,
        top_ok=(st <= transfer_tension_limit(fci)),
        bot_ok=(sb >= transfer_comp_limit(fci)),
    )


def batched_transfer(Pi: float, n_batch: int, n_total: int, sec: Section,
                     e: float, fci: float, M_sw_kNm: float = 0.0) -> StageStress:
    """支架上分批張拉：本批張 n_batch/n_total 組鋼腱，自重未活化（M_sw 預設 0）。

    過平衡頂緣引張隨張拉組數成正比；分批可把瞬時頂緣拉應力壓回容許值內，
    其餘鋼腱待脫架（自重活化）後再張。
    """
    return stage_stress(Pi * n_batch / n_total, sec, e, M_sw_kNm, fci)


# ── H3 平衡懸臂工法（Balanced Cantilever）──
# 對應算例_懸臂工法施工階段設計、公式卡_懸臂工法施工階段設計。

def variable_depth(x_from_mid: float, h_pier: float, h_mid: float,
                   half_span: float) -> float:
    """變深度箱梁拋物線斷面高 h(x) = h_mid + (h_pier − h_mid)·(x/半跨)²。

    x_from_mid：自跨中量起的距離（與 h_pier/h_mid/half_span 同單位，如 m）。
    端點：x=0 → h_mid；x=half_span → h_pier。
    """
    return h_mid + (h_pier - h_mid) * (x_from_mid / half_span) ** 2


def cantilever_moment(weights, arms, ft_load: float = 0.0, ft_arm: float = 0.0) -> float:
    """逐步懸臂彎矩 = Σ(G_i·arm_i) + 掛籃 G_FT·arm_FT（忠實對力臂求和）。

    weights/arms：各已澆節塊自重與其對驗核斷面的力臂（等長序列）；
    ft_load/ft_arm：掛籃（Form Traveler）重量與力臂。
    ⚠️ 力臂參考點由呼叫者決定，需自重項與掛籃項一致；算例公布的 94,025 其掛籃項
       以墩 CL 為準、自重項以 0 號塊端(x=4m)為準（混用，偏保守約 4%）——詳 README。
    """
    return sum(w * a for w, a in zip(weights, arms)) + ft_load * ft_arm


def long_term_deflection(d_elastic: float, phi: float) -> float:
    """長期下撓 ≈ δ_elastic·(1 + φ)（潛變放大；懸臂端典型 φ≈2.0）。"""
    return d_elastic * (1 + phi)
