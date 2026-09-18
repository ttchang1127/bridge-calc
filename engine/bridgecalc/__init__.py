"""bridgecalc — 橋梁設計計算引擎（單一真理源）。

對應 網頁計算器SOP §十一（混合架構）：本套件為 Python 權威計算核心，
驅動 handcalcs 計算書與 JS 網頁前端，並以知識庫算例 results 當回歸測試黃金答案。

階段 1（本版）：sections / prestress（含非線性損失耦合）/ loads / service / allowables。
"""
from .model import Section, Tendon
from .blister import (blister_local_bearing, blister_bursting, blister_tieback,
                      blister_spalling, blister_interface_shear,
                      blister_geometry_check, blister_design, BlisterDesign,
                      LocalBearingResult, BurstResult, TiebackResult,
                      SpallResult, InterfaceResult, BlisterGeomResult)
from .tendon_profile import (duct_size_check, DuctSizeResult, TW_DUCT_MAX_ID,
                             DUCT_AREA_RATIO)
from .staging import (redistribution_factor, creep_redistribution,
                      simple_span_dl_moment, span_by_span_dead_load,
                      redistribution_is_linear, prestress_M2_redistribution,
                      timing_sensitivity, simple_span_tendon_segs, staged_envelope,
                      simple_span_dl_shear, staged_shear_row,
                      StagedEnvRow, RedistFactor, RedistPoint,
                      StagingResult, M2RedistResult, TimingRow)
from .prestress import (compute_losses, LossResult, loss_profile,
                        LossProfileResult, LossPoint, parabolic_e, udl_moment)
from .loads import combinations, lane_live_load
from .service import stresses, Pe_min_zero_tension
from .shear import (shear_web, shear_web_at, stirrup_max_spacing_TW, stirrup_pick_spacing,
                    stirrup_zones, STD_STIRRUP_SPACINGS, ShearResult, phiVn, Av_s_min_TW,
                    principal_tension_limit_TW)
from .flexure import flexural_strength, FlexureResult, beta1
from .deflection import deflection_analysis, DeflectionResult
from .influence import (il_moment_simple, il_shear_simple, il_moment_peak,
                        max_moment_moving, abs_max_moment, lane_moment_simple,
                        hl93_per_lane_moment, moment_envelope_simple,
                        taiwan_per_lane_moment, taiwan_per_lane_shear, taiwan_impact,
                        taiwan_lane_moment, taiwan_lane_shear,
                        taiwan_truck_moment, taiwan_truck_shear, taiwan_moment_envelope)
from .fatigue import fatigue_check, FatigueResult, stirrup_fatigue
from .torsion import torsion_check, TorsionResult
from .transverse import slab_flexure, SlabResult, As_min_slab
from .temperature import (temp_gradient_AASHTO, ThermalBand, ThermalResult,
                          self_equilibrating_stress, thermal_service_check)
from .bearing import bearing_check, BearingResult
from .anchorage import anchorage_check, AnchorageResult, spiral_local_bearing
from .expansion import expansion_joint, JointResult
from .influence_cont import (cont_support_moments_point, cont_support_moments_uniform,
                             cont_moment_il, cont_dl_moment, cont_impact_length,
                             taiwan_lane_reduction, taiwan_cont_live_moment, ContLiveMoment,
                             taiwan_cont_envelope, ContEnvelopeRow,
                             cont_shear_il, cont_dl_shear, cont_shear_impact_length,
                             ContLiveShear, taiwan_cont_live_shear, ContShearRow,
                             taiwan_cont_shear_at, secondary_shear, design_shear_with_V2,
                             taiwan_rear_spacings, taiwan_cont_shear_envelope, ShearScanRow,
                             cont_shear_design_scan)
from .variable_section import (section_from_dims, HaunchProfile, haunch_profile,
                               ContFlex, gauss_integrate)
from .continuous import (secondary_moment, primary_moment, flexural_strength_T,
                        TFlexureResult, pier_service_stress,
                        ParabolaSeg, parabola_seg, ContTendonProfile, cont_tendon_segs,
                        TendonGroup, primary_moment_at, group_breaks,
                        ForceMethodM2Result, secondary_moments_force, continuous_prestress,
                        PierCapTendonResult, pier_cap_tendon_segs, groups_prestress_at)
from .tendon_profile import (tendon_profile, TendonProfileResult, equivalent_load,
                            end_slope, radius_of_curvature, balance_ratio, friction_loss,
                            duct_layout, duct_spacing_required, DuctLayoutResult,
                            top_slab_tendon_check, TopSlabTendonResult,
                            anchor_slip_loss, AnchorSlipResult, pier_cap_tendon_force, CapTendonForce,
                            tendon_slip_loss, segmented_tendon_force, SegmentedTendonForce,
                            segmented_friction_profile,
                            parabolic_curv_segs, friction_angle, friction_at,
                            friction_profile, FrictionProfileResult,
                            tendon_forces, assign_jack, TendonForceResult)
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
from .seismic import (min_falloff_length, ground_relative_displacement,
                     required_falloff_length, restrainer_yield_strength, EPSILON_G,
                     isolation_bilinear, damping_correction_B1, damping_correction_BS,
                     isolation_design, IsolationResult,
                     overstrength_moment, capacity_shear, rho_s_circular,
                     Ash_rectangular, confinement_spacing_limit, plastic_hinge_length,
                     liquefaction_reduction_DE)
from .retrofit import (cracked_na_depth, cracked_inertia, initial_concrete_strain,
                      cfrp_km1, cfrp_allowable_strain, xi_fb,
                      cfrp_moment_capacity, CFRPResult,
                      plate_moment_capacity, plate_dev_length, PlateResult,
                      enlargement_moment_capacity, EnlargeResult)
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
    "taiwan_truck_moment", "taiwan_truck_shear", "taiwan_moment_envelope",
    "fatigue_check", "FatigueResult", "stirrup_fatigue",
    "torsion_check", "TorsionResult",
    "slab_flexure", "SlabResult", "As_min_slab", "temp_gradient_AASHTO", "ThermalBand", "ThermalResult", "self_equilibrating_stress", "thermal_service_check", "bearing_check", "BearingResult", "anchorage_check", "AnchorageResult", "spiral_local_bearing", "expansion_joint", "JointResult",
    "secondary_moment", "primary_moment", "flexural_strength_T", "TFlexureResult", "pier_service_stress",
    "cont_support_moments_point", "cont_support_moments_uniform", "cont_moment_il",
    "cont_dl_moment", "cont_impact_length", "taiwan_lane_reduction", "taiwan_cont_live_moment",
    "ContLiveMoment", "taiwan_cont_envelope", "ContEnvelopeRow",
    "cont_shear_il", "cont_dl_shear", "cont_shear_impact_length", "ContLiveShear",
    "taiwan_cont_live_shear", "ContShearRow", "taiwan_cont_shear_at", "secondary_shear",
    "design_shear_with_V2", "shear_web_at", "taiwan_rear_spacings", "taiwan_cont_shear_envelope",
    "ShearScanRow", "cont_shear_design_scan", "stirrup_max_spacing_TW", "stirrup_pick_spacing",
    "stirrup_zones", "STD_STIRRUP_SPACINGS",
    "section_from_dims", "HaunchProfile", "haunch_profile", "ContFlex", "gauss_integrate",
    "ParabolaSeg", "parabola_seg", "ContTendonProfile", "cont_tendon_segs",
    "TendonGroup", "primary_moment_at", "group_breaks",
    "ForceMethodM2Result", "secondary_moments_force", "continuous_prestress",
    "PierCapTendonResult", "pier_cap_tendon_segs", "groups_prestress_at", "TopSlabTendonResult", "top_slab_tendon_check",
    "anchor_slip_loss", "AnchorSlipResult", "pier_cap_tendon_force", "CapTendonForce",
    "tendon_slip_loss", "segmented_tendon_force", "SegmentedTendonForce", "segmented_friction_profile",
    "tendon_profile", "TendonProfileResult", "equivalent_load", "end_slope",
    "radius_of_curvature", "balance_ratio", "friction_loss",
    "duct_layout", "duct_spacing_required", "DuctLayoutResult",
    "parabolic_curv_segs", "friction_angle", "friction_at",
    "friction_profile", "FrictionProfileResult",
    "tendon_forces", "assign_jack", "TendonForceResult",
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
    "min_tendon_groups", "required_drape", "min_section_modulus_Sb",
    "min_falloff_length", "ground_relative_displacement", "required_falloff_length",
    "restrainer_yield_strength", "EPSILON_G",
    "isolation_bilinear", "damping_correction_B1", "damping_correction_BS",
    "isolation_design", "IsolationResult",
    "overstrength_moment", "capacity_shear", "rho_s_circular", "Ash_rectangular",
    "confinement_spacing_limit", "plastic_hinge_length", "liquefaction_reduction_DE",
    "cracked_na_depth", "cracked_inertia", "initial_concrete_strain",
    "cfrp_km1", "cfrp_allowable_strain", "xi_fb", "cfrp_moment_capacity", "CFRPResult",
    "plate_moment_capacity", "plate_dev_length", "PlateResult",
    "enlargement_moment_capacity", "EnlargeResult",
    "blister_local_bearing", "blister_bursting", "blister_tieback",
    "blister_spalling", "blister_interface_shear", "blister_geometry_check",
    "blister_design", "BlisterDesign", "LocalBearingResult", "BurstResult",
    "TiebackResult", "SpallResult", "InterfaceResult", "BlisterGeomResult",
    "redistribution_factor", "creep_redistribution", "simple_span_dl_moment",
    "span_by_span_dead_load", "redistribution_is_linear",
    "prestress_M2_redistribution", "timing_sensitivity",
    "simple_span_tendon_segs", "staged_envelope", "StagedEnvRow",
    "simple_span_dl_shear", "staged_shear_row",
    "RedistFactor", "RedistPoint", "StagingResult", "M2RedistResult", "TimingRow",
    "duct_size_check", "DuctSizeResult", "TW_DUCT_MAX_ID", "DUCT_AREA_RATIO",
    "allowables",
]
