"""服務性應力（C1）與最小預力反解。壓為負（compression negative）。

對應公式卡 C1（服務性應力限制）、算例_後張箱梁服務性應力驗算。
"""
from .model import Section


def stresses(Pe: float, section: Section, e: float, M_kNm: float):
    """跨中頂/底緣應力（MPa，壓為負）。Pe 為 N，M_kNm 為 kN·m。"""
    M = M_kNm * 1e6
    sb = -Pe / section.A - Pe * e / section.Sb + M / section.Sb
    st = -Pe / section.A + Pe * e / section.St - M / section.St
    return st, sb


def Pe_min_zero_tension(section: Section, e: float, M_kNm: float) -> float:
    """★ 符號反解（驗算變設計）：底緣零拉所需最小有效預力 Pe（N）。

    由 σ_bot = −Pe/A − Pe·e/Sb + M/Sb ≤ 0 解得：
        Pe,min = M·A / (Sb + A·e)
    （SymPy 可符號導出此閉合式，見 Python計算引擎_handcalcs範例.py）
    """
    M = M_kNm * 1e6
    return M * section.A / (section.Sb + section.A * e)
