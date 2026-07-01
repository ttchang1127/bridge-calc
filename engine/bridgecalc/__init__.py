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
                        hl93_per_lane_moment, moment_envelope_simple,
                        taiwan_per_lane_moment, taiwan_per_lane_shear, taiwan_impact,
                        taiwan_lane_moment, taiwan_lane_shear,
                        taiwan_truck_moment, taiwan_moment_envelope)
from .fatigue import fatigue_check, FatigueResult, stirrup_fatigue
from .torsion import torsion_check, TorsionResult
from .transverse import slab_flexure, SlabResult, As_min_slab
from .temperature import (temp_gradient_AASHTO, ThermalBand, ThermalResult,
                          self_equilibrating_stress, thermal_service_check)
from .bearing import bearing_check, BearingResult
from .anchorage import anchorage_check, AnchorageResult, spiral_local_bearing
from .expansion import expansion_joint, JointResult
from .continuous import (secondary_moment, primary_moment, flexural_strength_T,
                        TFlexureResult, pier_service_stress)
from .tendon_profile import (tendon_profile, TendonProfileResult, equivalent_load,
                            end_slope, radius_of_curvature, balance_ratio, friction_loss)
from .stm import (general_zone_burst, STMResult, burst_force, burst_depth, f_cu,
                 strut_capacity, node_capacity, tie_reinforcement,
                 BETA_NODE, BETA_STRUT)
from .durability import (grout_qc_check, GroutQCResult, rebar_stress_limit,
                        rebar_stress_ok, pc_fatigue_limit, design_life,
                        GROUT, REBAR_LIMIT, DESIGN_LIFE)
from .construction import (stage_stress, batched_transfer, StageStress,
                          transfer_tension_limit, transfer_comp_limit,
                          variable_depth, cantilever_moment, long_term_deflection)
from .launching import (launching_cantilever_moment, launching_span_moment,
                       centric_prestress_required, launching_bottom_stress,
                       n_tendons, jacking_force, bearing_stress)
from .segmental import (segment_weight, joint_min_prestress, joint_compression,
                       shear_key_design_capacity, shear_key_utilization, bonded_pt_ratio,
                       JOINT_MIN_COMPRESSION_MPa, BONDED_PT_MIN_RATIO)
from .design import min_tendon_groups, required_drape, min_section_modulus_Sb
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
    "hl93_per_lane_moment", "moment_envelope_simple",
    "taiwan_per_lane_moment", "taiwan_per_lane_shear", "taiwan_impact",
    "taiwan_lane_moment", "taiwan_lane_shear",
    "taiwan_truck_moment", "taiwan_moment_envelope",
    "fatigue_check", "FatigueResult", "stirrup_fatigue",
    "torsion_check", "TorsionResult",
    "slab_flexure", "SlabResult", "As_min_slab", "temp_gradient_AASHTO", "ThermalBand", "ThermalResult", "self_equilibrating_stress", "thermal_service_check", "bearing_check", "BearingResult", "anchorage_check", "AnchorageResult", "spiral_local_bearing", "expansion_joint", "JointResult",
    "secondary_moment", "primary_moment", "flexural_strength_T", "TFlexureResult", "pier_service_stress",
    "tendon_profile", "TendonProfileResult", "equivalent_load", "end_slope",
    "radius_of_curvature", "balance_ratio", "friction_loss",
    "general_zone_burst", "STMResult", "burst_force", "burst_depth", "f_cu",
    "strut_capacity", "node_capacity", "tie_reinforcement", "BETA_NODE", "BETA_STRUT",
    "grout_qc_check", "GroutQCResult", "rebar_stress_limit", "rebar_stress_ok",
    "pc_fatigue_limit", "design_life", "GROUT", "REBAR_LIMIT", "DESIGN_LIFE",
    "stage_stress", "batched_transfer", "StageStress",
    "transfer_tension_limit", "transfer_comp_limit",
    "variable_depth", "cantilever_moment", "long_term_deflection",
    "launching_cantilever_moment", "launching_span_moment", "centric_prestress_required",
    "launching_bottom_stress", "n_tendons", "jacking_force", "bearing_stress",
    "segment_weight", "joint_min_prestress", "joint_compression",
    "shear_key_design_capacity", "shear_key_utilization", "bonded_pt_ratio",
    "JOINT_MIN_COMPRESSION_MPa", "BONDED_PT_MIN_RATIO",
    "min_tendon_groups", "required_drape", "min_section_modulus_Sb", "allowables",
]
