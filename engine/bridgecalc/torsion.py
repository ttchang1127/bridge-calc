"""扭力（D2）：開裂扭矩 T_cr 與忽略門檻。

對應公式卡 D2、算例_扭力設計。箱梁閉合斷面 T_cr 極大，偏載扭矩常 < 門檻 → 可免顯式
扭矩設計，但仍須配置閉合箍筋。單位：力 N、長度 mm、應力 MPa；T 以 kN·m 表示。
"""
from dataclasses import dataclass
from math import sqrt
from .model import Section


@dataclass
class TorsionResult:
    fpc: float        # 形心預力壓應力 MPa
    Tcr: float        # 開裂扭矩 kN·m
    threshold: float  # 忽略門檻 0.25·φ·Tcr kN·m
    neglect: bool     # Tu < 門檻 → 可忽略顯式扭矩設計
    need_closed_stirrup: bool  # 箱梁恆需閉合箍筋


def torsion_check(section: Section, Pe: float, fc: float,
                  Acp: float, pc: float, Tu_kNm: float, phi: float = 0.90) -> TorsionResult:
    """AASHTO 開裂扭矩與門檻判斷。

    Acp：外周封閉面積 mm²；pc：外周長 mm（箱梁取毛斷面外輪廓）。
    """
    fpc = Pe / section.A
    k = 0.125 * sqrt(fc)                              # N/mm²
    Tcr = k * (Acp ** 2 / pc) * sqrt(1 + fpc / k) / 1e6   # kN·m
    threshold = 0.25 * phi * Tcr
    return TorsionResult(fpc, Tcr, threshold,
                         Tu_kNm < threshold, True)
