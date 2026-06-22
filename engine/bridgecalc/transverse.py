"""橫向設計（D3）：頂板／懸臂翼板 RC 板撓曲（每公尺寬）。

對應公式卡 D3、算例_橫向設計、橋面板設計法快查卡。
單筋矩形 RC 撓曲：a=As·fy/(0.85f'c·b)，φMn=φ·As·fy·(d−a/2)。
亦適用一般橋面板撓曲檢核。單位：mm / MPa；M 以 kN·m/m 表示（b=1000）。
"""
from dataclasses import dataclass


@dataclass
class SlabResult:
    a: float       # 等效矩形深度 mm
    phiMn: float   # 設計彎矩強度 kN·m/m
    ok: bool       # φMn ≥ Mu？


def slab_flexure(Mu_kNm: float, As: float, d: float, fc: float, fy: float,
                 b: float = 1000.0, phi: float = 0.9) -> SlabResult:
    """每公尺寬 RC 板撓曲設計強度。As：mm²/m；d：有效深 mm。"""
    a = As * fy / (0.85 * fc * b)
    phiMn = phi * As * fy * (d - a / 2) / 1e6
    return SlabResult(a, phiMn, phiMn >= Mu_kNm)


def As_min_slab(fc: float, fy: float, b: float, d: float) -> float:
    """最小鋼筋 0.03·(f'c/fy)·b·d（mm²/m）。"""
    return 0.03 * fc / fy * b * d
