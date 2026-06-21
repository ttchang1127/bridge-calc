"""極限強度彎曲（M1）：AASHTO 近似法 fps、中性軸 c、Mn、延性 φ、Mcr 下限。

對應公式卡 M1（極限強度彎曲設計）、算例_極限強度彎曲設計。
僅含矩形壓力區（NA 在翼板內）；c > hf 時設 in_flange=False 提醒需 T 斷面公式。
單位：力 N、長度 mm、應力 MPa；Mn/Mu/Mcr 以 kN·m 表示。
"""
from dataclasses import dataclass
from math import sqrt
from .model import Section, Tendon


@dataclass
class FlexureResult:
    beta1: float
    c: float          # 中性軸深度 mm
    in_flange: bool   # NA 在翼板內？（否則需 T 斷面）
    a: float          # 等效矩形深度 mm
    fps: float        # 鋼腱極限應力 MPa
    Mn: float         # 公稱彎矩強度 kN·m
    eps_t: float      # 最外受拉應變
    phi: float        # 強度折減係數
    phiMn: float      # 設計彎矩強度 kN·m
    CR: float         # 強度比 φMn/Mu
    Mcr: float        # 開裂彎矩 kN·m
    lower_bound: float  # min(1.33Mu, 1.2Mcr)
    ok: bool


def beta1(fc: float) -> float:
    """AASHTO 等效矩形應力區係數。"""
    if fc <= 28:
        return 0.85
    return max(0.65, 0.85 - 0.05 * (fc - 28) / 7)


def flexural_strength(tendon: Tendon, section: Section, fc: float,
                      b_eff: float, hf: float, dp: float,
                      Mu_kNm: float, Pe: float, e: float,
                      dt: float = None) -> FlexureResult:
    """跨中正彎矩極限強度（矩形壓力區，後張無黏結/有黏結近似法）。"""
    Aps = tendon.Aps
    fpu = tendon.fpu
    fpy = 0.90 * fpu                       # 低鬆弛
    k = 2 * (1.04 - fpy / fpu)             # = 0.28
    b1 = beta1(fc)

    c = Aps * fpu / (0.85 * fc * b1 * b_eff + k * Aps * fpu / dp)
    in_flange = c <= hf
    a = b1 * c
    fps = fpu * (1 - k * c / dp)
    Mn = Aps * fps * (dp - a / 2) / 1e6    # kN·m

    dt = dt if dt is not None else dp
    eps_t = (dt - c) / c * 0.003
    if eps_t >= 0.005:
        phi = 1.0                          # 拉力控制（PC）
    elif eps_t <= 0.002:
        phi = 0.75
    else:
        phi = 0.75 + 0.25 * (eps_t - 0.002) / 0.003   # 過渡區
    phiMn = phi * Mn
    CR = phiMn / Mu_kNm

    fr = 0.97 * sqrt(fc)                                # 破裂模數（本庫採值）
    fcpe = Pe / section.A + Pe * e / section.Sb         # 底緣預壓 MPa
    Mcr = section.Sb * (fr + fcpe) / 1e6               # kN·m
    lower = min(1.33 * Mu_kNm, 1.2 * Mcr)

    return FlexureResult(b1, c, in_flange, a, fps, Mn, eps_t, phi, phiMn, CR,
                         Mcr, lower, phiMn >= Mu_kNm and phiMn >= lower)
