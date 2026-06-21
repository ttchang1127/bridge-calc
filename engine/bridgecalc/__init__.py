"""bridgecalc — 橋梁設計計算引擎（單一真理源）。

對應 網頁計算器SOP §十一（混合架構）：本套件為 Python 權威計算核心，
驅動 handcalcs 計算書與 JS 網頁前端，並以知識庫算例 results 當回歸測試黃金答案。

階段 1（本版）：sections / prestress（含非線性損失耦合）/ loads / service / allowables。
"""
from .model import Section, Tendon
from .prestress import compute_losses, LossResult
from .loads import combinations, lane_live_load
from .service import stresses, Pe_min_zero_tension
from .shear import (shear_web, ShearResult, phiVn, Av_s_min_TW,
                    principal_tension_limit_TW)
from .flexure import flexural_strength, FlexureResult, beta1
from .deflection import deflection_analysis, DeflectionResult
from .influence import (il_moment_simple, il_shear_simple, il_moment_peak,
                        max_moment_moving, abs_max_moment, lane_moment_simple,
                        hl93_per_lane_moment, moment_envelope_simple)
from . import allowables

__all__ = [
    "Section", "Tendon", "compute_losses", "LossResult",
    "combinations", "lane_live_load", "stresses", "Pe_min_zero_tension",
    "shear_web", "ShearResult", "phiVn", "Av_s_min_TW",
    "principal_tension_limit_TW",
    "flexural_strength", "FlexureResult", "beta1",
    "deflection_analysis", "DeflectionResult",
    "il_moment_simple", "il_shear_simple", "il_moment_peak",
    "max_moment_moving", "abs_max_moment", "lane_moment_simple",
    "hl93_per_lane_moment", "moment_envelope_simple", "allowables",
]
