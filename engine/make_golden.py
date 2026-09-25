"""重新產生 golden_answers.json（Python 引擎與 JS 前端的共用驗證源）。

CI 在每次 push 時執行此檔，確保 golden 永遠來自最新的 bridgecalc。
執行：python3 make_golden.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bridgecalc import (Section, Tendon, box_slab_thickness_TW, compute_losses, combinations,
                        lane_live_load, stresses, Pe_min_zero_tension,
                        shear_web, phiVn, flexural_strength, deflection_analysis,
                        il_moment_peak, abs_max_moment, lane_moment_simple,
                        hl93_per_lane_moment, moment_envelope_simple,
                        taiwan_per_lane_moment, taiwan_per_lane_shear, taiwan_truck_shear, taiwan_impact,
                        taiwan_truck_moment, taiwan_lane_moment, taiwan_moment_envelope,
                        fatigue_check, stirrup_fatigue, torsion_check,
                        slab_flexure, As_min_slab, temp_gradient_AASHTO,
                        bearing_check, anchorage_check, spiral_local_bearing, expansion_joint,
                        ThermalBand, self_equilibrating_stress, thermal_service_check,
                        secondary_moment, primary_moment, flexural_strength_T, pier_service_stress,
                        tendon_profile, general_zone_burst, node_capacity, f_cu,
                        blister_stm, node_capacity_aashto, strut_fcu_aashto,
                        grout_qc_check, rebar_stress_limit, pc_fatigue_limit, design_life, GROUT,
                        batched_transfer, stage_stress, transfer_tension_limit,
                        variable_depth, cantilever_moment, long_term_deflection,
                        launching_cantilever_moment, launching_span_moment,
                        centric_prestress_required, launching_bottom_stress,
                        n_tendons, jacking_force, bearing_stress,
                        segment_weight, joint_min_prestress, joint_compression,
                        shear_key_design_capacity, shear_key_utilization, bonded_pt_ratio,
                        min_tendon_groups, required_drape, min_section_modulus_Sb,
                        duct_layout, duct_spacing_required,
                        parabolic_curv_segs, friction_angle, friction_at, friction_profile,
                        tendon_forces, assign_jack,
                        parabola_seg, cont_tendon_segs, TendonGroup, primary_moment_at,
                        continuous_prestress, taiwan_cont_envelope, taiwan_cont_live_moment,
                        cont_moment_il, cont_dl_moment, taiwan_lane_reduction,
                        pier_cap_tendon_segs, top_slab_tendon_check,
                        section_from_dims, haunch_profile, ContFlex,
                        cont_shear_il, cont_dl_shear, taiwan_cont_live_shear, taiwan_cont_shear_at,
                        secondary_shear, design_shear_with_V2, shear_web_at, taiwan_rear_spacings,
                        anchor_slip_loss, pier_cap_tendon_force, cont_shear_design_scan,
                        groups_prestress_at, stirrup_max_spacing_TW, tendon_slip_loss,
                        segmented_tendon_force, segmented_friction_profile,
                        loss_profile, parabolic_e, udl_moment, blister_design,
                        redistribution_factor, creep_redistribution,
                        span_by_span_dead_load, redistribution_is_linear,
                        prestress_M2_redistribution, timing_sensitivity,
                        simple_span_tendon_segs, staged_envelope, duct_size_check,
                        positive_moment_connection, durability_cover,
                        aashto_creep, staging_phi, box_volume_surface, timing_sensitivity_aashto)
from bridgecalc.tendon_profile import duct_layout_bundled, end_zone_duct_check, web_duct_check, deviation_force_check, adjacent_duct_radial_check
from bridgecalc.staging import (pos_conn_strand, pos_conn_rebar, strand_stress_extended,
                                strand_length_required, negative_moment_connection,
                                neg_top_tension_limit, multi_stage_redistribution,
                                StageSpec, span_by_span_schedule, stage_moments_span_by_span,
                                cantilever_units, cantilever_layout, cantilever_M_I, cantilever_M_I_points,
                                cantilever_unbalanced, cantilever_falsework, cantilever_X1X0,
                                cantilever_schedule, cantilever_stages)
from bridgecalc import allowables
from bridgecalc import seismic as seis
from bridgecalc import retrofit as retro

sec = Section(5.065e6, 3.287e12, 1329, 2100)
ten = Tendon(8, 19, 1109)                       # 台灣 HS20-44 最小設計（HL-93 才需增配 21 股）


def _cont_case_groups(e_top=646.0):
    """算例_連續梁次彎矩：底板腱全長（各跨拋物線 端/墩 −80、跨中 +950）＋頂板腱 x=25~55（錨 +300、墩 −e_top）。"""
    cb = -(950 + 80) / 20 ** 2
    bot = TendonGroup(23700, [parabola_seg(0, 40, 20, 950, cb), parabola_seg(40, 80, 60, 950, cb)])
    ct = (300 + e_top) / 15 ** 2
    top = TendonGroup(12557, [parabola_seg(25, 40, 40, -e_top, ct), parabola_seg(40, 55, 40, -e_top, ct)])
    return bot, top


def _cont_P_e(groups, x):
    """斷面 x 的總預力 P 與合力偏心（形心下為正）；錨碇點以該處存在的鋼腱計。"""
    P = Pe = 0.0
    for g in groups:
        for sg in g.segs:
            if sg.x1 - 1e-9 <= x <= sg.x2 + 1e-9:
                P += g.P
                Pe += g.P * sg.e(x)
                break
    return P, Pe / P


def _cont_service_scan(groups, rows, fm, with_M2=True):
    """跨 1 逐斷面（含墩）服務性：回傳 (底緣最拉值, x, 頂緣最拉值, x)；M_ext 取正負包絡、加 M2。"""
    wb = wt = None
    for r in rows:
        if r.x > 40 + 1e-9:
            break
        P, e = _cont_P_e(groups, r.x)
        M2 = fm.M2_at(r.x) if with_M2 else 0.0
        for Ms in (r.Ms_pos, r.Ms_neg):
            st, sb = pier_service_stress(P * 1e3, sec, e, Ms + M2)
            if wb is None or sb > wb[0]:
                wb = (sb, r.x)
            if wt is None or st > wt[0]:
                wt = (st, r.x)
    return wb, wt


def _continuous_pier():
    bot, top = _cont_case_groups()
    g = [bot, top]
    r = continuous_prestress([40, 40], g)
    M2p = r.X[1]
    w_dc = sec.A / 1e6 * 24.5
    rows = taiwan_cont_envelope([40, 40], w_dc, 20, 2, n_per_span=40)
    pier = next(rw for rw in rows if abs(rw.x - 40) < 1e-9)
    e_eff = -(23700 * 80 + 12557 * 646) / 36257
    st_p, sb_p = zip(*(pier_service_stress(36257e3, sec, e_eff, Ms + M2p) for Ms in (pier.Ms_pos, pier.Ms_neg)))
    sb_p0 = pier_service_stress(36257e3, sec, e_eff, pier.Ms_neg)[1]
    (wb, xb), (wt, xt) = _cont_service_scan(g, rows, r)
    (wb0, xb0), _ = _cont_service_scan(g, rows, r, with_M2=False)
    Mu = -(pier.Mu_neg + M2p)
    Mu0 = -pier.Mu_neg
    ft = flexural_strength_T(4 * 2660, 1860, 40, 1400, 200, 700, 1975, Mu)
    ft0 = flexural_strength_T(4 * 2660, 1860, 40, 1400, 200, 700, 1975, Mu0)
    pos = max((rw for rw in rows if rw.x <= 40), key=lambda rw: rw.Mu_pos + r.M2_at(rw.x))
    tp = cont_tendon_segs([40, 40], -80, 950, -600, 0.12, 0.40)
    ga = [TendonGroup(23700, tp.segs), top]
    ra = continuous_prestress([40, 40], ga)
    (wba, xba), _ = _cont_service_scan(ga, rows, ra)
    return {
        "M1_mid_kNm": round(primary_moment_at(g, 20)), "M1_pier_kNm": round(primary_moment_at(g, 40)),
        "M2_x15_kNm": round(r.M2_at(15)), "M2_mid_kNm": round(r.M2_at(20)), "M2_pier_kNm": round(M2p),
        "M2_pier_bot_only_kNm": round(continuous_prestress([40, 40], [bot]).X[1]),
        "M2_pier_top_only_kNm": round(continuous_prestress([40, 40], [top]).X[1]),
        "int_M1_m_kNm2": round(r.b[0]), "flex_F_m": round(r.F[0][0], 3),
        "w_dc_kNm": round(w_dc, 2),
        "pier_M_dc_kNm": round(pier.M_dc), "pier_M_dw_kNm": round(pier.M_dw), "pier_M_ll_neg_kNm": round(pier.M_ll_neg),
        "pier_impact": round(pier.I_neg, 4),
        "e_pier_eff_mm": round(e_eff, 1),
        "pier_service_sigma_top_MPa": round(max(st_p), 2), "pier_service_sigma_bot_MPa": round(min(sb_p), 2),
        "pier_service_sigma_bot_noM2_MPa": round(sb_p0, 2),
        "span_max_sigma_bot_MPa": round(wb, 2), "span_max_sigma_bot_x_m": xb,
        "span_max_sigma_top_MPa": round(wt, 2), "span_max_sigma_top_x_m": xt,
        "span_max_sigma_bot_noM2_MPa": round(wb0, 2), "span_max_sigma_bot_noM2_x_m": xb0,
        "pos_Mu_max_kNm": round(pos.Mu_pos + r.M2_at(pos.x)), "pos_Mu_max_x_m": pos.x,
        "pier_Mu_kNm": round(Mu), "pier_Mu_noM2_kNm": round(Mu0),
        "pier_c_mm": round(ft.c), "pier_flanged": ft.flanged,
        "pier_fps_MPa": round(ft.fps), "pier_Mn_kNm": round(ft.Mn), "pier_phiMn_kNm": round(ft.phiMn),
        "pier_CR": round(ft.CR, 2), "pier_CR_noM2": round(ft0.CR, 2), "pier_inadequate": not ft.ok,
        "adj_M2_pier_kNm": round(ra.X[1]), "adj_span_max_sigma_bot_MPa": round(wba, 2),
        "_note": "2026-09-17 二次校正：外力改引擎實算（taiwan_cont_envelope：自重 A×24.5、SDL 20、HS20-44 2車道、"
                 "負彎矩車道 2 集中載重、衝擊依跨長）。算例原估活載 12,000/−15,000 → 實算 5,510/−5,966。"
                 "M2 力法全長正彎矩；全線服務性通過；中墩 Mu 含 M2 → CR 1.05✓（不計 M2 為 0.67✗）",
    }


def _pier_cap_tendon():
    """雙系統配束：墩頂局部腱線形＋頂板構造檢核＋與全長腱合成的 M2。"""
    pc = pier_cap_tendon_segs([40, 40], 15, 300, -646)
    top = TendonGroup(12557, pc.segs)
    chk75 = top_slab_tendon_check(2100, 1329, 250, -646, 4, 5800 - 700, 100, 75)
    chk40 = top_slab_tendon_check(2100, 1329, 250, -646, 4, 5800 - 700, 100, 40)
    bad = top_slab_tendon_check(2100, 1329, 250, -700, 12, 5800 - 700, 100, 40, 25, "od")
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    dual = continuous_prestress([40, 40], [TendonGroup(23724, tp.segs), top])
    pc3 = pier_cap_tendon_segs([30, 40, 30], 20, 300, -600)
    return {
        "c_mm_per_m2": round(pc.segs[0].c, 6), "R_min_mm": round(pc.R_min), "anchors_m": pc.anchors,
        "M2_top_only_kNm": round(continuous_prestress([40, 40], [top]).X[1]),
        "slab75_e_hi": chk75.e_hi, "slab75_e_lo": chk75.e_lo, "slab75_ok": chk75.ok,
        "slab40_e_hi": chk40.e_hi, "slab40_e_lo": chk40.e_lo, "slab40_s_clear": chk40.s_clear,
        "bad_in_slab": bad.in_slab, "bad_top_cover": bad.top_cover, "bad_s_clear": round(bad.s_clear, 2), "bad_s_ok": bad.s_ok,
        "dual_M2_pier_kNm": round(dual.X[1], 1),
        "three_span_lengths": pc3.lengths,
        "_note": "算例_連續梁次彎矩頂板腱（b 15、錨 +300、墩 −646）；頂板 250＝75+100+75 恰容一層",
    }


def _variable_section_haunch():
    """變斷面（中墩底板加厚）：斷面性質、變 EI 柔度（階梯 EI 解析驗證）、M2 與包絡。"""
    dims = (11000, 250, 5800, 200, 350, 2, 2100)
    ref = section_from_dims(*dims)
    pr = haunch_profile([40, 40], *dims, 400, 8)
    pr0 = haunch_profile([40, 40], *dims, 200, 8)
    s40 = pr.section_at(40)
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    g = [TendonGroup(23724, tp.segs)]
    bot, top = _cont_case_groups()
    e_h = taiwan_cont_envelope([40, 40], 124.0925, 20, 2, n_per_span=20, stiff=pr)
    e_0 = taiwan_cont_envelope([40, 40], 124.0925, 20, 2, n_per_span=20, stiff=pr0)
    step = ContFlex([40, 40], lambda x: 2.0 if abs(x - 40) < 8 else 1.0, breaks=[32, 48])
    return {
        "ref_A": round(ref.A), "ref_I_e12": round(ref.I / 1e12, 5), "ref_yb": round(ref.yb, 3),
        "pier_A": round(s40.A), "pier_I_e12": round(s40.I / 1e12, 5), "pier_yb": round(s40.yb, 3),
        "dyb_pier": round(pr.dyb(40), 3), "I_rel_pier": round(pr.I_rel(40), 6), "bot_t_x36": pr.bot_t_at(36),
        "numeric0_M2_pier": round(continuous_prestress([40, 40], g, stiff=pr0).X[1], 2),
        "numeric0_env_pier_Mu_neg": round(e_0[20].Mu_neg, 1),
        "step_EI_dl_pier_w10": round(step.support_moments_dist(lambda t: 10.0)[1], 3),
        "step_EI_point_p20": round(step.support_moments_point(20)[1], 6),
        "haunch_M2_pier_default": round(continuous_prestress([40, 40], g, stiff=pr).X[1], 1),
        "haunch_M2_pier_dual": round(continuous_prestress([40, 40], [bot, top], stiff=pr).X[1], 1),
        "haunch_env_pier_M_dc": round(e_h[20].M_dc, 1), "haunch_env_pier_M_ll_neg": round(e_h[20].M_ll_neg, 1),
        "haunch_env_pier_Mu_neg": round(e_h[20].Mu_neg, 1), "haunch_env_x16_Mu_pos": round(e_h[8].Mu_pos, 1),
        "_note": "底板 200→400（墩心）線性加厚、單側 8 m；階梯 EI（墩兩側 8 m 內 2 倍）解析值 −2,406.349（w=10）",
    }


def _cont_shear_taiwan():
    """連續梁剪力：影響線、DL、HS20-44 活載、次剪力 V2 與算例墩左側 d_v 斷面腹板抗剪。"""
    sp = [40, 40]
    bot, top = _cont_case_groups()
    fm = continuous_prestress(sp, [bot, top])
    x = 40 - 1.512                                   # d_v = max(0.9 d_p, 0.72h) = 0.72×2100（d_p=ȳb−e 較小）
    P = Pe = Ps = 0.0
    for g in (bot, top):
        for sg in g.segs:
            if sg.x1 <= x <= sg.x2:
                P += g.P
                Pe += g.P * sg.e(x)
                Ps += g.P * sg.slope(x) / 1000
                break
    row = taiwan_cont_shear_at(sp, x, "L", sec.A / 1e6 * 24.5, 20, 2)
    V2 = secondary_shear(fm, sp, x, "L")
    Vu = design_shear_with_V2(row.Vu_pos, row.Vu_neg, V2)
    sh = shear_web_at(P * 1e3, sec, Ps / P, 40, 250, 1512, Vu / 2 * 1e3, 2)
    lvL = taiwan_cont_live_shear(sp, 40, "L")
    lvS = taiwan_cont_live_shear([40], 0, "R")
    lv3 = taiwan_cont_live_shear([30, 40, 30], 30, "R")
    return {
        "il_x10R_p30": round(cont_shear_il(sp, 10, 30, "R"), 8), "il_pierL_p20": round(cont_shear_il(sp, 40, 20, "L"), 8),
        "il_pierR_p60": round(cont_shear_il(sp, 40, 60, "R"), 8),
        "dl_pierL_w10": round(cont_dl_shear(sp, 10, 40, "L"), 6), "dl_end_w10": round(cont_dl_shear(sp, 10, 0, "R"), 6),
        "pierL_truck_neg": round(lvL.truck_neg, 2), "pierL_lane_neg": round(lvL.lane_neg, 2),
        "single_truck": round(lvS.truck_pos, 2), "single_lane": round(lvS.lane_pos, 2),
        "three_pier1R_lane_pos": round(lv3.lane_pos, 2), "three_pier1R_I": round(lv3.I, 6),
        "case_x_m": x, "case_V_dc": round(row.V_dc, 2), "case_V_ll_neg": round(row.V_ll_neg, 2),
        "case_Vu_neg": round(row.Vu_neg, 2), "case_V2": round(V2, 2), "case_Vu_design": round(Vu, 2),
        "case_P_slope_kN": round(Ps, 2), "case_Vp_per_web_kN": round(sh.Vp / 1e3, 2),
        "case_sigma1": round(sh.sigma1, 3), "case_Vcw_kN": round(sh.Vcw / 1e3, 2), "case_Av_s_req": round(sh.Av_s_req, 4),
        "case_phiVn_D16x2_s250_kN": round(phiVn(sh.Vcw, 397.4 / 250, 1512) / 1e3, 1),
        "case_phiVn_D16x2_s200_kN": round(phiVn(sh.Vcw, 397.4 / 200, 1512) / 1e3, 1),
        "_note": "算例雙系統 40+40 墩左 d_v 斷面：每腹板 Vu 2,876 kN；D16×2@250 φVn 2,766✗、@200 2,980✓；V2 +434 於墩左有利→不折減",
    }


def _rear_axle_scan():
    """HS20-44 中後軸距 4.25～9.15 m 掃描：短跨連續梁墩頂由長軸距控制；簡支跨恆為 4.25。"""
    from bridgecalc.influence_cont import _ILCache, _truck_extremes, _GRID
    c = _ILCache([15, 15])
    N = int(round(c.tot / _GRID))
    vals = [c.eta(15, k * _GRID) for k in range(N + 1)]
    fixed = _truck_extremes(vals, c.tot, (15, lambda p, sd: c.eta(15, p)), spacings=[4.25])
    lv = taiwan_cont_live_moment([15, 15], 15)
    s40 = taiwan_cont_live_moment([40], 20)
    return {
        "spacings_n": len(taiwan_rear_spacings()), "spacings_last": taiwan_rear_spacings()[-1],
        "p15_truck_neg_scan": round(lv.truck_neg, 3), "p15_V_neg": lv.V_neg,
        "p15_truck_neg_fixed425": round(fixed[1], 3),
        "simple40_x20_truck": round(s40.truck_pos, 3), "simple40_x20_V": s40.V_pos,
        "_note": "15+15 墩頂：兩後軸分落兩跨負區 → V=9.15 控制；簡支 40 m 仍 4.25（與 influence.py 相同）",
    }


def _simple_shear_dv():
    """簡支 40 m 近支承 d_v 斷面設計剪力（analyzer ⑤ 自動帶入 Vu）：DL 解析＋HS20-44（衝擊長度＝至較遠支點）。"""
    L, x = 40.0, 1.692
    w_dc = sec.A / 1e6 * 24.5
    r = taiwan_cont_shear_at([L], x, "R", w_dc, 20, 2)
    lv = taiwan_cont_live_shear([L], x, "R")
    return {
        "x_m": x, "V_dc": round(r.V_dc, 3), "V_dw": round(r.V_dw, 3), "truck": round(lv.truck_pos, 3),
        "lane": round(lv.lane_pos, 3), "I": round(r.I, 6), "V_ll": round(r.V_ll_pos, 3),
        "Vu_total": round(r.Vu_pos, 2), "Vu_per_web": round(r.Vu_pos / 2, 2),
        "_note": "analyzer ⑤ Vu 自動帶入值；算例_腹板抗剪 2,329 係等值均布活載 42.5 kN/m＋SDL 因數 1.25 之近似",
    }


def _pier_cap_tendon_losses():
    """頂板腱逐點損失（摩擦＋錨具滑移 6 mm＋ES＋長期）與對中墩 M2／強度的影響（算例 §3.3 對照）。"""
    pc = pier_cap_tendon_segs([40, 40], 15, 300, -646)
    Aps = 4 * 2660.0
    fcgp = 36257e3 / sec.A + 36257e3 * (-276.0) * (-646.0) / sec.I - (-24819e6) * (-646.0) / sec.I
    ES = 3 / 8 * 7.33 * fcgp
    other = ES + 155.0
    f40, f25, f385 = (pier_cap_tendon_force(pc, x, 1395.0, Aps, other) for x in (40, 25, 38.488))
    bot, _ = _cont_case_groups()
    top = TendonGroup(lambda x: pier_cap_tendon_force(pc, x, 1395.0, Aps, other).P, pc.segs)
    r = continuous_prestress([40, 40], [bot, top])
    rows = taiwan_cont_envelope([40, 40], sec.A / 1e6 * 24.5, 20, 2, n_per_span=40)
    pier = next(rw for rw in rows if abs(rw.x - 40) < 1e-9)
    Mu = -(pier.Mu_neg + r.X[1])
    ft = flexural_strength_T(Aps, 1860, 40, 1400, 200, 700, 1975, Mu)
    sl = anchor_slip_loss(0, 15, 5.0, 6)
    return {
        "fcgp_top_pier": round(fcgp, 3), "ES": round(ES, 2), "other": round(other, 2),
        "pier_friction": round(f40.friction, 2), "pier_slip": round(f40.slip, 2), "pier_fpe": round(f40.fpe, 2),
        "pier_P_kN": round(f40.P, 1), "anchor_slip": round(f25.slip, 2), "anchor_fpe": round(f25.fpe, 2),
        "dv_x_slip": round(f385.slip, 3), "dv_x_P_kN": round(f385.P, 1),
        "M2_pier_kNm": round(r.X[1], 1), "pier_Mu_kNm": round(Mu, 1), "pier_CR": round(ft.CR, 3),
        "slip_capped_dsigma": round(sl.dsigma, 3), "slip_capped": sl.capped,
        "_note": "算例假設 f_pe 1,180（短期 60 估）；精算墩頂 f_pe≈1,126（摩擦 103＋ES 12＋長期 155）→ P 12,558→約 11,975；"
                 "M2 上升（頂板腱 M2 為負）、A_ps 改實際 10,640 → CR 約 1.04",
    }


def _cont_shear_scan():
    """算例雙系統 40+40 全長剪力掃描＋箍筋分區（D16×2 肢；錨碇點前後加取樣）。"""
    bot, top = _cont_case_groups()
    fm = continuous_prestress([40, 40], [bot, top])
    rows, zones = cont_shear_design_scan([40, 40], groups_prestress_at([bot, top]), fm, sec.A / 1e6 * 24.5, 20, 2,
                                         2100, 1329, 40, 250, 2, 397.4, step=1.0,
                                         extra=[24.95, 25.05, 54.95, 55.05], sec=sec)
    by = {round(r.x, 2): r for r in rows}
    a_l, a_r = by[24.95], by[25.05]
    worst = max(rows, key=lambda r: r.Av_s_req)
    return {
        "n_rows": len(rows), "worst_x": round(worst.x, 3), "worst_Av_s": round(worst.Av_s_req, 4),
        "worst_s": worst.s_pick, "worst_halved": worst.halved,
        "anchorL_Vp": round(a_l.Vp_web, 1), "anchorR_Vp": round(a_r.Vp_web, 1),
        "anchorL_s": a_l.s_pick, "anchorR_s": a_r.s_pick,
        "zones": [[round(z[0], 3), round(z[1], 3), z[2]] for z in zones],
        "max_spacing_halved": stirrup_max_spacing_TW(1e6, 40, 250, 1512, 2100)[0],
        "max_spacing_normal": stirrup_max_spacing_TW(5e5, 40, 250, 1512, 2100)[0],
        "_note": "中墩兩側須 @200（V_s 超 0.33√f'c·b·d_v → 上限減半 300）；頂板腱錨碇左側 V_p 驟減 → 需 @400，右側 @500",
    }


def _tendon_slip_G1():
    """全長腱錨具滑移（Δa 6 mm）：簡支 40 m 雙端張拉跨中不受影響（驗證 compute_losses「跨中＝0」）；
    連續 80 m 交錯張拉端部 d_v 斷面 P 下降。"""
    fpj = 0.75 * 1860
    seg40 = parabolic_curv_segs(40, 1109)
    p40 = fpj * friction_at(20, 40, seg40, 0.25, 0.003, "start") / 20
    Ls40 = (6 * 195000 / 1000 / p40) ** 0.5
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    cs = [(g.x1, g.x2, 2 * abs(g.c) / 1000) for g in tp.segs]
    lay = duct_layout(8, 2, sec.yb - 1109, y_b=sec.yb, h=2100)
    jk = assign_jack(lay.n_per_web, "alt")
    ts, d = [], 0
    for _iw in range(2):
        d = 0
        for (y, cnt) in lay.rows:
            for j in range(cnt):
                if d < lay.n_per_web:
                    ts.append({"y": y, "x": lay.x_offsets[j], "jack": jk[d]})
                    d += 1
    lo = compute_losses(Tendon(8, 19, 1109), sec, 0, 0)
    oth = lo.total - lo.friction
    f0 = tendon_forces(ts, 1.5, 80, cs, sec.yb, fpj, 19 * 140, oth)
    f6 = tendon_forces(ts, 1.5, 80, cs, sec.yb, fpj, 19 * 140, oth, slip_mm=6)
    m0 = tendon_forces(ts, 40, 80, cs, sec.yb, fpj, 19 * 140, oth, slip_mm=6)
    return {
        "simple40_p_MPa_per_m": round(p40, 4), "simple40_L_set_m": round(Ls40, 3),
        "simple40_mid_unaffected": Ls40 < 20,
        "simple40_slip_anchor": round(tendon_slip_loss(0, 40, seg40, fpj, "both", slip_mm=6), 2),
        "simple40_slip_mid": round(tendon_slip_loss(20, 40, seg40, fpj, "both", slip_mm=6), 2),
        "cont80_slip_start_jack_x1_5": round(tendon_slip_loss(1.5, 80, cs, fpj, "start", slip_mm=6), 2),
        "cont80_slip_end_jack_x1_5": round(tendon_slip_loss(1.5, 80, cs, fpj, "end", slip_mm=6), 2),
        "cont80_Pe_x1_5_noslip_kN": round(f0.Pe_total / 1e3, 1), "cont80_Pe_x1_5_slip6_kN": round(f6.Pe_total / 1e3, 1),
        "cont80_drop_pct": round((f6.Pe_total / f0.Pe_total - 1) * 100, 2),
        "cont80_Pe_pier_slip6_kN": round(m0.Pe_total / 1e3, 1),
        "_note": "Δa 6 mm。簡支 40 m 雙端：L_set 14.1 < 半長 20 → 跨中確實為 0（compute_losses 假設成立）；"
                 "連續 80 m 交錯：端部 d_v 斷面 P 降 9.1%（半數腱在該端受滿滑移 153 MPa）→ V_p 同比下降",
    }


def _mid_anchor_G1():
    """中間錨碇段（G1 不合格決策樹 >60 m）：40+40 全長腱 80 m 摩擦超 20% → 墩頂設中間錨碇。"""
    fpj = 0.75 * 1860
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    cs = [(g.x1, g.x2, 2 * abs(g.c) / 1000) for g in tp.segs]
    one0 = segmented_friction_profile([0, 80], cs, fpj, 0.25, 0.003, "both", 0)
    one6 = segmented_friction_profile([0, 80], cs, fpj, 0.25, 0.003, "both", 6)
    two6 = segmented_friction_profile([0, 40, 80], cs, fpj, 0.25, 0.003, "both", 6)
    two0 = segmented_friction_profile([0, 40, 80], cs, fpj, 0.25, 0.003, "both", 0)
    three6 = segmented_friction_profile([0, 80 / 3, 160 / 3, 80], cs, fpj, 0.25, 0.003, "both", 6)
    f_mid = segmented_tendon_force(40, [0, 40, 80], cs, fpj, 8 * 19 * 140, 0.0, 0.25, 0.003, "both", 6)
    return {
        "one_seg_max_pct": round(one0["max_ratio"] * 100, 2), "one_seg_x_max": round(one0["x_max"], 1),
        "one_seg_slip6_max_pct": round(one6["max_ratio"] * 100, 2),
        "two_seg_slip6_max_pct": round(two6["max_ratio"] * 100, 2), "two_seg_slip6_x_max": round(two6["x_max"], 1),
        "two_seg_fric_only_max_pct": round(two0["max_ratio"] * 100, 2),
        "three_seg_slip6_max_pct": round(three6["max_ratio"] * 100, 2),
        "pier_fpe_two_seg": round(f_mid.fpe, 1), "pier_friction_two_seg": round(f_mid.friction, 1),
        "pier_slip_two_seg": round(f_mid.slip, 1),
        "_note": "一段 80 m 雙端張拉：摩擦 20.28% @墩頂（超 G1 20% 門檻，雙端張拉救不了）；"
                 "墩頂設中間錨碇（2×40 m）→ 摩擦＋滑移 14.16%（純摩擦 10.36%）；"
                 "再分三段 13.92% 幾無改善——每增一個錨碇就多一組滑移損失",
    }


def _cont_envelope_taiwan():
    """連續梁解析影響線＋台灣 HS20-44 包絡的檢核點。"""
    sp = [40, 40]
    lv_p = taiwan_cont_live_moment(sp, 40)
    lv_15 = taiwan_cont_live_moment(sp, 15)
    lv_s = taiwan_cont_live_moment([40], 20)
    sp3 = [30, 40, 30]
    lv3 = taiwan_cont_live_moment(sp3, 30)
    rows = taiwan_cont_envelope(sp, 124.0925, 20, 2, n_per_span=20)
    return {
        "il_pier_p20": round(cont_moment_il(sp, 40, 20), 6),
        "il_x15_p15": round(cont_moment_il(sp, 15, 15), 6),
        "dl_pier_w130_4": round(cont_dl_moment(sp, 130.4, 40), 3), "dl_x15_w130_4": round(cont_dl_moment(sp, 130.4, 15), 3),
        "pier_truck_neg": round(lv_p.truck_neg, 2), "pier_lane_neg": round(lv_p.lane_neg, 2), "pier_I_neg": round(lv_p.I_neg, 6),
        "x15_truck_pos": round(lv_15.truck_pos, 2), "x15_lane_pos": round(lv_15.lane_pos, 2),
        "x15_truck_neg": round(lv_15.truck_neg, 2), "x15_lane_neg": round(lv_15.lane_neg, 2),
        "single_span_lane_mid": round(lv_s.lane_pos, 2), "single_span_truck_mid": round(lv_s.truck_pos, 2),
        "three_span_pier1_lane_neg": round(lv3.lane_neg, 2), "three_span_pier1_truck_neg": round(lv3.truck_neg, 2),
        "three_span_pier1_I_neg": round(lv3.I_neg, 6),
        "env_min_Mu_neg_kNm": round(min(r.Mu_neg for r in rows), 1),
        "env_max_Mu_pos_kNm": round(max(r.Mu_pos for r in rows), 1),
        "env_max_Ms_pos_x_m": max(rows, key=lambda r: r.Ms_pos).x,
        "lane_reduction_3": taiwan_lane_reduction(3), "lane_reduction_4": taiwan_lane_reduction(4),
        "_note": "三彎矩解析（EI常數）。§3.9 負彎矩車道 2 集中載重（不同跨最負縱距）；§3.13 衝擊 L：正＝該跨、負＝相鄰兩跨平均",
    }


def _cont_tendon_force_default():
    """analyzer ⑦ 預設全長連續腱（40+40、e_end 0、e_mid 1109、e_pier −600、k 0.12、r 0.42、P 常數 23,724）。"""
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    r = continuous_prestress([40, 40], [TendonGroup(23724, tp.segs)])
    tp3 = cont_tendon_segs([30, 40, 30], 0, 900, -500)
    r3 = continuous_prestress([30, 40, 30], [TendonGroup(20000, tp3.segs)])
    return {
        "n_segs": len(tp.segs), "infl_m": [round(v, 3) for v in tp.infl], "lows_m": [round(v, 3) for v in tp.lows],
        "R_min_mm": round(tp.R_min), "e_at_10": round(tp.e_at(10), 2), "e_at_36": round(tp.e_at(36), 2),
        "M2_pier_kNm": round(r.X[1], 1), "M2_x20_kNm": round(r.M2_at(20), 1),
        "M1_pier_kNm": round(primary_moment_at([TendonGroup(23724, tp.segs)], 40), 1),
        "three_span_X_kNm": [round(v, 1) for v in r3.X],
        "_note": "force method (EI const). analyzer IL 等效載重法 μ=K=0 時 13,937@P23,307 → ×23,724/23,307 = 14,186 互驗",
    }
M_DC, M_DW = 24800, 4000
M_LL = lane_live_load(taiwan_per_lane_moment(40), 2, 1.0)   # HS20-44 2 車道 = 6,837
L = compute_losses(ten, sec, M_DC, M_DW)
c = combinations(M_DC, M_DW, M_LL)
st, sb = stresses(L.Pe, sec, ten.e, c["Service_I"])
# d_v 斷面設計剪力：引擎實算（DL 解析＋HS20-44 取大、衝擊長度＝至較遠支點）÷2 腹板 ≈ 2,298
# （原為「按 HS20/HL-93 縮放 ≈ 2,069」；算例_腹板抗剪的 2,329 係等值均布活載 42.5 kN/m＋SDL 因數 1.25 之近似）
Vu_HS20 = taiwan_cont_shear_at([40], 1.692, "R", sec.A / 1e6 * 24.5, 20, 2).Vu_pos / 2
sh = shear_web(L.Pe, sec, ten.e, 40, 250, 1692, Vu_HS20 * 1e3, 1692, 40000)
fx = flexural_strength(ten, sec, 40, 8000, 250, 1880, c["Strength_I"], L.Pe, ten.e)
w_LL_HS20 = 56.7 * taiwan_per_lane_moment(40) / hl93_per_lane_moment(40)   # 撓度等效 UDL 按比例 ≈ 34.2
df = deflection_analysis(40000, 29700, sec, 144, L.Pe, ten.e, w_LL_HS20)
R_LL_HS20 = 290 * taiwan_per_lane_shear(40) / 588          # 支承活載反力 HS20/HL-93 縮放 ≈ 179
an = anchorage_check(ten.Pi / 1e3, 8, 260, 2100, 4)
sp = spiral_local_bearing(an.Pu, 2919, 8.47, 104044, 50, 380)   # 螺旋圍束 D16@50

# ── 線形/STM/耐久/施工 擴充（本批，對齊各 verified 算例）──
ten21 = Tendon(8, 21, 1109)                                     # G1 HL-93 敘述軌
w_DL_g1 = 8 * (M_DC + M_DW) * 1e6 / 40000 ** 2                  # = 144 kN/m
g1 = tendon_profile(ten21.Pi, compute_losses(ten21, sec, M_DC, M_DW).Pe, 1109, 40000, w_DL_g1)
# G1 管道實配：參考橋 8 組 / 2 腹板 / φ100 / 腹板 350
# dl  ＝台灣規範明文（§8.25.1 保護層 40、§8.25.2 淨距 max(40,1.5d_agg)）→ 兩列兩層，e=1109 排得下
# dls ＝保守側（PTBG 實務保護層 75、淨距再取 ≥孔道外徑）→ 單列四層，同一配置就排不下
dl = duct_layout(8, 2, sec.yb - 1109, duct_od=100, web_t=350, cover=40, s_v=40, d_agg=25,
                 y_b=sec.yb, h=2100, rule="tw")
dls = duct_layout(8, 2, sec.yb - 1109, duct_od=100, web_t=350, cover=75, s_v=100, d_agg=25,
                  y_b=sec.yb, h=2100, rule="od")
# 摩擦損失沿長度分佈：四種張拉端配置（參考橋 L=40 m、a=1109）
fsg = parabolic_curv_segs(40, 1109)
fps = {jk: friction_profile(40, fsg, jack=jk) for jk in ('both', 'start', 'end', 'alt')}
# 逐腱合成：兩種排法 × 交錯端，看合力位置是否偏離幾何形心（x=L/4 損失左右差最大處）
def _tf(rows, mode):
    jk = assign_jack(len(rows), mode)
    ts = [{"no": "T%d" % (i + 1), "y": rows[i][0], "x": rows[i][1], "jack": jk[i]}
          for i in range(len(rows))]
    return tendon_forces(ts, 10, 40, fsg, sec.yb, 1395, 140 * 19, 180.0)
_rowsA = [(-80, 0), (120, 0), (320, 0), (520, 0)]        # 單列四層：序號＝層號
_rowsB = [(150, -85), (150, 85), (290, -85), (290, 85)]  # 兩列兩層：同層左右配對
tfA, tfB = _tf(_rowsA, "alt"), _tf(_rowsB, "alt")
tfS = _tf(_rowsA, "start")
dlp = duct_layout(8, 2, sec.yb + 900, duct_od=100, web_t=350, cover=40, s_v=40, d_agg=25,
                  y_b=sec.yb, h=2100, rule="tw")      # 墩頂 e_pier = −900（形心上方）
f2 = general_zone_burst(14850e3, 1050, 2100, 150, L.Pe, sec.A, 540000, 12830e3, 40, 420)
h8 = batched_transfer(29700e3, 8, 8, sec, 1109, 32)            # S2 全 PT
h4 = batched_transfer(29700e3, 4, 8, sec, 1109, 32)            # S2 分批 4 組
h3s = stage_stress(29700e3, sec, 1109, 24800, 32)             # S3 脫架
w_h3 = [643, 599, 550, 497, 439, 385, 353, 341]
a_h3 = [x - 4 for x in [6.75, 11.25, 15.75, 20.25, 24.75, 29.25, 33.75, 38.25]]
Mpos_h4 = launching_span_moment(120.0, 40.0)
Pc_h4 = centric_prestress_required(Mpos_h4, 3.093e9, 4.870e6, 1.5)
Vkey_h6 = shear_key_design_capacity(350, 0.439)
# 反解設計庫（階段 4）：設計=驗算之逆，閉環對齊參考橋
Pe_min_ref = Pe_min_zero_tension(sec, ten.e, c["Service_I"])
d_ngroups = min_tendon_groups(Pe_min_ref, 19, L.fpe)
d_drape = required_drape(0.985, w_DL_g1, 40000, compute_losses(ten21, sec, M_DC, M_DW).Pe)
d_Sbmin = min_section_modulus_Sb(L.Pe, sec.A, ten.e, c["Service_I"], 0.0)

# 耐震閉式（台灣單軌）：S2 隔震迭代結果算一次
seis_S2 = seis.isolation_design(8000, 400, 6000, 0.60)

# 補強閉式（中國 JTG 單軌）：R1/R2/R4 同梁 b400/h800/h01750/As1964/C30
_a = 6.667                                          # α_Es = E_s/E_c = 2.0e5/3.0e4
retro_x1 = retro.cracked_na_depth(400, 750, 1964, _a)
retro_Icr = retro.cracked_inertia(400, retro_x1, 1964, 750, _a)
retro_ec1 = retro.initial_concrete_strain(200, retro_x1, 3.0e4, retro_Icr)
retro_epsf = retro.cfrp_allowable_strain(2, 2.4e5, 0.167, 0.0155)
retro_R1 = retro.cfrp_moment_capacity(400, 800, 750, 13.8, 0.0033, 1964, 330,
                                      100.2, 2.4e5, retro_epsf, 0.0003)
retro_R2 = retro.plate_moment_capacity(400, 800, 750, 13.8, 0.0033, 1964, 330,
                                       800, 305, 2.06e5, retro_x1, retro_ec1)
retro_R4 = retro.enlargement_moment_capacity(400, 860, (1964*750+982*860)/2946,
                                             13.8, 0.0033, 1964, 330, 982, 330,
                                             2.0e5, retro_x1, retro_ec1)


def _loss_profile_G1():
    """沿長度的預力 P(x)（取代單點 fric_ratio）：40 m 參考橋、Δa 6 mm 雙端張拉。

    鎖住三件單點式看不到的事：①端部摩擦 0 但滑移最大 ②滑移在 L_set 外歸零
    ③最不利斷面隨是否計滑移而搬家（跨中→端部）。
    """
    L, a = 40.0, 870.0
    segs = parabolic_curv_segs(L, a)
    MDf = udl_moment(sec.A / 1e6 * 24.5, L)
    efn = parabolic_e(ten.e, ten.e - a, L)
    lp0 = loss_profile(ten, sec, L, efn, MDf, segs, mu=0.25, K=0.003,
                       jack="both", slip_mm=0.0, n=40)
    lp1 = loss_profile(ten, sec, L, efn, MDf, segs, mu=0.25, K=0.003,
                       jack="both", slip_mm=6.0, n=40)
    p0, pm = lp1.at(0.0), lp1.at(20.0)
    single = compute_losses(ten, sec, MDf(L / 2), 0.0,
                            fric_ratio=friction_at(L / 2, L, segs, 0.25, 0.003, "both"))
    return {
        "config": "40m參考橋 e_mid1109/垂度870/μ0.25/K0.003/雙端張拉/Δa6mm/n=40",
        "end_fric_pct": round(p0.fric / ten.fpj * 100, 2),
        "end_slip_pct": round(p0.slip / ten.fpj * 100, 2),
        "end_loss_pct": round(p0.loss_pct * 100, 2),
        "end_Pe_kN": round(p0.Pe / 1e3, 1),
        "mid_fric_pct": round(pm.fric / ten.fpj * 100, 2),
        "mid_slip_pct": round(pm.slip / ten.fpj * 100, 2),
        "mid_loss_pct": round(pm.loss_pct * 100, 2),
        "mid_Pe_kN": round(pm.Pe / 1e3, 1),
        "x_Pemin_noslip_m": round(lp0.x_Pemin, 1),
        "x_Pemin_slip_m": round(lp1.x_Pemin, 1),
        "Pe_min_noslip_kN": round(lp0.Pe_min / 1e3, 1),
        "Pe_min_slip_kN": round(lp1.Pe_min / 1e3, 1),
        "Pe_avg_slip_kN": round(lp1.Pe_avg / 1e3, 1),
        "single_point_mid_Pe_kN": round(single.Pe / 1e3, 1),
        "single_point_overestimate_at_1m_pct":
            round((single.Pe - lp1.at(1.0).Pe) / lp1.at(1.0).Pe * 100, 2),
        "_note": "單點式(compute_losses取跨中)全長都用 23,687kN；實際端部僅 23,565kN(滑移11.47%但摩擦0)。"
                 "無滑移時 Pe_min 在跨中(摩擦控制)、計滑移後搬到端部——最不利斷面會搬家是單點式看不到的。"
                 "跨中滑移=0 印證 compute_losses 註解的假設(L_set 14.1m < 半長 20m)。"}


def _blister_B1():
    """中間錨碇齒塊錨碇區（算例_外置PT補強_錨固齒塊中間錨碇設計）：12φ15.2、P_s 1,764 kN。

    齒塊與端部錨碇／轉向塊的根本差異：後方無鋼腱延伸提供回力 → Tie-back 不得省略；
    介面剪力逼近全腱力（cos5°）而非偏向力（sin5°），相差約 11 倍。
    """
    d = blister_design(Ps_kN=1764, Pu_kN=3000, Pd_kN=3000,
                       a_plate=200, b_plate=200, L_b=400, W_b=300, D_b=250,
                       A_bearing=60000, fci=35, f_cb=8.0, A_cb=40000,
                       alpha_deg=5, mu=1.0, fy=420, fsd=360, fc=40, straight_have=400)
    return {
        "config": "12φ15.2/P_s1764/P_u3000/P_d3000/錨板200x200/齒塊400x300x250/f'ci35/f_cb8",
        "f_b_MPa": round(d.bearing.f_b, 1),
        "f_b_allow_MPa": round(d.bearing.f_b_allow, 2),
        "bearing_ok": d.bearing.ok,
        "spiral_factor_req": round(d.bearing.spiral_factor_req, 2),
        "F_burst_kN": round(d.burst.F_burst, 1),
        "d_burst_mm": round(d.burst.d_burst, 1),
        "As_burst_mm2": round(d.burst.As_burst, 1),
        "T_tieback_kN": round(d.tie.T_req, 1),
        "C_precomp_kN": round(d.tie.C_precomp, 1),
        "fs_allow_MPa": round(d.tie.fs_allow, 1),
        "As_tie_mm2": round(d.tie.As, 1),
        "As_tie_conservative_mm2": round(d.tie.As_conservative, 1),
        "F_spall_kN": round(d.spall.F_spall, 2),
        "As_spall_mm2": round(d.spall.As_spall, 1),
        "V_interface_kN": round(d.face.V_int, 1),
        "As_vf_mm2": round(d.face.As_vf, 1),
        "tau_interface_MPa": round(d.face.tau, 3),
        "tau_cap_MPa": round(d.face.tau_cap, 2),
        "V_cap_kN": round(d.face.V_cap, 1),
        "interface_area_ok": d.face.area_ok,
        "A_interface_req_mm2": round(d.face.A_req, 0),
        "governing": d.governing,
        "As_total_conservative_mm2": round(d.As_total_conservative, 1),
        "geom_ok": d.geom.ok,
        "_note": "局部承壓 75.0 > 25.5 不通過 → 必須配錨具螺旋筋(需提升2.94倍)，這是齒塊錨板面積受限的必然結果不是錯誤。"
                 "配筋由介面剪力摩擦控制(4,881mm²)——齒塊傳 P·cos5°=1,757kN，轉向塊只傳 P·sin5°=154kN，差11倍。"
                 "Tie-back 兩值並列：含既有預壓抵扣 488、不抵扣 1,778；補強既有橋時預壓分布難確認應取後者，引擎不替使用者選邊。"
                 "🔴算例漏檢介面面積上限：τ=17.57 > min(0.25f'c,10.3)=10.0 MPa，400×250 的介面差1.76倍，"
                 "加多少鋼筋都沒用只能加大齒塊(600×300 即通過，且配筋需求完全不變)。K1/K2 隨規範版次不同。"}


def _staging_redist_S1():
    """施工階段體系轉換的潛變重分配：40+40 逐跨施工（先簡支後連續）。

    逐跨施工的橋既不是簡支也不是連續——墩頂拿到七成連續負彎矩，跨中卻仍保留
    超過連續值的正彎矩，**兩頭都要設計**。這是實務常見的控制因素。
    """
    sp, w, dphi = [40.0, 40.0], 5.065 * 24.5, 1.7
    xs, m1, m2 = span_by_span_dead_load(sp, w, 40)
    r = creep_redistribution(xs, m1, m2, dphi, 0.8, "trost")
    rd = creep_redistribution(xs, m1, m2, dphi, 0.8, "dischinger")
    # 簡支時張拉的腱只能是各跨自己的拋物線（連續腱線形只在合龍後才存在）
    cp = continuous_prestress(sp, [TendonGroup(P=23725.0, segs=simple_span_tendon_segs(sp, 1109.0))])
    xp = [float(i) for i in range(81)]
    rm = prestress_M2_redistribution(sp, xp, [cp.M2_at(x) for x in xp], dphi)
    tim = timing_sensitivity(2.0, [("7d", 0.25), ("28d", 0.55), ("180d", 1.35)],
                             0.0, -w * 1600 / 8)
    mid, pier = r.at(20.0), r.at(40.0)
    return {
        "config": "40+40 逐跨施工/w=124.09kN/m/Δφ=1.7/Trost χ=0.8",
        "lam_trost": round(r.factor.lam, 4),
        "lam_dischinger": round(rd.factor.lam, 4),
        "lam_trost_invalid_at_dphi10": not redistribution_factor(10.0, 0.8, "trost").valid,
        "mid_M_I_kNm": round(mid.M_I, 1),
        "mid_M_II_kNm": round(mid.M_II, 1),
        "mid_M_inf_kNm": round(mid.M_inf, 1),
        "mid_ratio_to_cont": round(mid.M_inf / mid.M_II, 3),
        "pier_M_I_kNm": round(pier.M_I, 1),
        "pier_M_II_kNm": round(pier.M_II, 1),
        "pier_M_inf_kNm": round(pier.M_inf, 1),
        "pier_ratio_to_cont": round(pier.M_inf / pier.M_II, 3),
        "restraint_linear": redistribution_is_linear(r, sp),
        "M2_pier_cont_kNm": round(rm.M2_pier_cont, 1),
        "M2_pier_inf_kNm": round(rm.M2_pier_inf, 1),
        "pier_total_cont_kNm": round(pier.M_II + rm.M2_pier_cont, 1),
        "pier_total_sbs_kNm": round(pier.M_inf + rm.M2_pier_inf, 1),
        "M2_pier_closed_form_Pa_kNm": round(23725.0 * 1109.0 / 1000, 1),
        "lam_7d": round(tim[0].lam, 4),
        "lam_180d": round(tim[2].lam, 4),
        "M_pier_7d_kNm": round(tim[0].M_pier, 1),
        "M_pier_180d_kNm": round(tim[2].M_pier, 1),
        "_note": "M(∞)=M_I+λ(M_II−M_I)。跨中 15,880＝連續值的128%(保留部分簡支正彎矩)、"
                 "墩頂 −17,878＝連續值的72%——逐跨施工把不利處從墩頂搬到跨中，兩頭都要設計。"
                 "🔑 束制彎矩 ΔM=λ(M_II−M_I) 不對應外載故支承間必為線性，這是本項最強的自我驗證。"
                 "M₂ 在簡支張拉時＝0(靜定)，連續後由潛變生成 λ·M₂,cont；🔴M₂,cont 須以各跨簡支拋物線算"
                 "(力法 26,311＝閉合解 P·a)，2026-09-18 首版誤用連續腱線形得 14,185。長期墩頂合計 +1,075 為**正彎矩**"
                 "→需正彎矩接頭(AASHTO 5.14.1.4.9a；C5.14.1.4.5 早合龍時梁持續上拱致連續橫隔梁底開裂)。"
                 "🔴 決定 λ 的是**剩餘**潛變 Δφ 而非 φ(∞)：7天合龍 λ=0.729 vs 180天 λ=0.428，"
                 "趕工早合龍反而使墩頂更不利(−18,097 vs −10,613，差71%)。"}


def _staged_envelope_S2():
    """體系轉換套進連續梁包絡（40+40、台灣活載 2 車道、DW 20 kN/m、Δφ 1.7 Trost）。

    兩狀態取不利：t₁（合龍當下，自重＝簡支）與 ∞（長期重分配）；恆載有利時取 γ_min。
    """
    sp, w, lam = [40.0, 40.0], 5.065 * 24.5, 1.7 / (1 + 0.8 * 1.7)
    rows = taiwan_cont_envelope(sp, w, 20.0, 2, 20)
    xs = [r.x for r in rows]
    ip = min(range(len(xs)), key=lambda i: abs(xs[i] - 40))
    cp = continuous_prestress(sp, [TendonGroup(P=23725.0, segs=simple_span_tendon_segs(sp, 1109.0))])
    M2 = [cp.M2_at(x) for x in xs]
    a = staged_envelope(sp, rows, w, lam)                              # 連續後張拉（M₂ 另計）
    b = staged_envelope(sp, rows, w, lam, M2=M2, ps_at_simple=True)    # 簡支時張拉
    im0 = max(range(len(rows)), key=lambda i: rows[i].Mu_pos)
    ima = max(range(len(a)), key=lambda i: a[i].Mu_pos)
    return {
        "config": "40+40/w124.09/DW20/2車道/Δφ1.7 Trost/γ_DC 1.25|0.90、γ_DW 1.50|0.65",
        "mono_mid_Mu_pos_kNm": round(rows[im0].Mu_pos, 1),
        "mono_pier_Mu_neg_kNm": round(rows[ip].Mu_neg, 1),
        "sbs_mid_Mu_pos_kNm": round(a[ima].Mu_pos, 1),
        "sbs_mid_x_m": round(a[ima].x, 2),
        "sbs_mid_gov": a[ima].gov_pos,
        "sbs_pier_Mu_neg_kNm": round(a[ip].Mu_neg, 1),
        "sbs_pier_gov_neg": a[ip].gov_neg,
        "pss_pier_Mu_neg_kNm": round(b[ip].Mu_neg, 1),
        "pss_pier_Mu_pos_kNm": round(b[ip].Mu_pos, 1),
        "pss_pier_gov_pos": b[ip].gov_pos,
        "_note": "逐跨施工使跨中 Mu⁺ 較一次成形大 43%，由 t₁ 控制(剛合龍時自重全為簡支值、活載已作用於連續體系)。"
                 "簡支時張拉：長期中墩 Mu⁺ +263＞0→需正彎矩接頭(AASHTO 5.14.1.4.9a：max(係數化正束制彎矩,0.6M_cr))；"
                 "自重對正彎矩有利須取 γ_min 0.90，若仍乘 1.25 得 −9,394 會誤判無需接頭。"
                 "AASHTO 5.14.1.4.2 束制彎矩有利時不得計入→兩狀態取不利即自動滿足。"}


def _duct_size_T83():
    """套管尺寸：台灣 §8.25.4（面積 ≥ 2 倍）、表 8.3（最大內徑）；PTI Table 4.3（2.25/2.5）。"""
    a = duct_size_check(19, 100.0)
    b = duct_size_check(19, 90.0, rule="pti_pull")
    c = duct_size_check(19, 105.0)
    return {
        "ratio_19x152_id100": round(a.ratio, 4), "area_ok_tw": a.area_ok,
        "id_max_tw_19x152_mm": a.id_max_tw, "od_gt_id_when_od_unknown": a.od_gt_id,
        "ratio_id90": round(b.ratio, 4), "area_ok_id90_pti_pull": b.area_ok,
        "id_ok_id105": c.id_ok,
        "_note": "表 8.3 與 PTI Table 4.4 列的都是**內徑**，外徑依廠商(PTI §4.4.5)；排列用外徑、面積比用內徑。"
                 "參考橋 φ100 是內徑，排列檢核把它當外徑用→偏不保守，須輸入廠商外徑。"}


def _cont_envelope_gamma_min():
    """連續梁包絡恆載因數依有利／不利取 max／min（AASHTO γ_p）：反曲點附近的反轉範圍。"""
    rows = taiwan_cont_envelope([40.0, 40.0], 5.065 * 24.5, 20.0, 2, 20)
    at = lambda x: min(rows, key=lambda r: abs(r.x - x))
    return {"x28_Mu_neg_kNm": round(at(28.0).Mu_neg, 1), "x32_Mu_pos_kNm": round(at(32.0).Mu_pos, 1),
            "pier_Mu_pos_kNm": round(at(40.0).Mu_pos, 1), "x20_Mu_neg_kNm": round(at(20.0).Mu_neg, 1),
            "_note": "2026-09-19 前恆載一律 1.25/1.50：x=28 Mu⁻ +1,529（判無負彎矩）→ −163；x=32 Mu⁺ −1,591 → +343；"
                     "墩頂 Mu⁺ −37,023 → −24,937。控制斷面（跨中正、墩頂負）不變，golden 既有項零變動。"}


def _staged_shear_S3():
    """體系轉換的恆載剪力＋箍筋分區（同 cont_shear_scan 雙系統案，Δφ 1.7 Trost）。"""
    lam = 1.7 / (1 + 0.8 * 1.7)
    bot, top = _cont_case_groups()
    fm = continuous_prestress([40, 40], [bot, top])
    args = ([40, 40], groups_prestress_at([bot, top]), fm, sec.A / 1e6 * 24.5, 20, 2,
            2100, 1329, 40, 250, 2, 397.4)
    kw = dict(step=1.0, extra=[24.95, 25.05, 54.95, 55.05], sec=sec)
    r0, z0 = cont_shear_design_scan(*args, **kw)
    r1, z1 = cont_shear_design_scan(*args, stage_lam=lam, **kw)
    return {
        "end_dv_Vu_web_mono": round(r0[0].Vu_web, 1), "end_dv_Vu_web_sbs": round(r1[0].Vu_web, 1),
        "end_zone_s_mono": z0[0][2], "end_zone_s_sbs": z1[0][2], "end_zone_to_sbs": round(z1[0][1], 3),
        "zones_sbs": [[round(z[0], 3), round(z[1], 3), z[2]] for z in z1],
        "_note": "逐跨施工剛合龍時自重剪力為簡支 wL/2(>連續 3wL/8)，端支承 d_v Vu/腹板 2,053→2,441(+19%)，"
                 "端部 0–4.5 m 箍筋 @450→@250。束制剪力 λ(V_II−V_I) 每跨常數(束制彎矩線性)。"}


def _pos_moment_conn_S4():
    """連續橫隔梁正彎矩接頭（AASHTO 5.14.1.4.9a／.4）：40+40 簡支張拉中墩 Mu⁺ 263（staged_envelope_S2）。"""
    a = positive_moment_connection(262.9, sec.I, sec.yb, 40.0, d=sec.h - 80)
    b = positive_moment_connection(262.9, sec.I, sec.yb, 40.0, age_days=90, d=sec.h - 80)
    return {"fr_MPa": round(a.fr, 4), "Mcr_kNm": round(a.Mcr, 1), "gov": a.governs, "M_req_kNm": round(a.M_req, 1),
            "As_est_mm2": round(a.As_est, 0), "gov_90d": b.governs, "M_req_90d_kNm": round(b.M_req, 1),
            "_note": "M_cr＝f_r·I_g/y_t：總毛斷面、橫隔梁混凝土 f_r 0.63√f'c、不計預力(C5.14.1.4.9a)。"
                     "實際正束制 263 ≪ 0.6M_cr 5,913 → 最小量控制；90 天簡化要 1.2M_cr 11,826——省計算不省鋼筋。"}


def _durability_cover_T12():
    """台灣第十二章耐久性保護層（表 12.2／12.5，p.341–343），取 max(表列, §8.25.1 基本 40 mm)。"""
    g50 = durability_cover("general", 50, "II", 0.45)
    g100 = durability_cover("general", 100, "II", 0.45)
    s50 = durability_cover("salt", 50, member="梁腹版外露面", zone="嚴重", wc=0.40)
    s100 = durability_cover("salt", 100, member="梁腹版外露面", zone="極嚴重", wc=0.40)
    return {"gen50_II_wc045_table": g50.table, "gen50_II_wc045_req": g50.required,
            "gen100_II_wc045_req": g100.required, "gen50_I_wc045_table": durability_cover("general", 50, "I", 0.45).table,
            "salt50_web_severe_req": s50.required, "salt100_web_extreme_req": s100.required,
            "deck_top_salt50_moderate": durability_cover("salt", 50, member="橋面版頂層筋", zone="中度", wc=0.45).required,
            "_note": "橋面板頂面 50 mm 只在鹽害中度 50 年成立（另 §7.1.5 表 7.2 RC 值）；一般環境外露面 50 年 w/c 0.45 表列 35 < §8.25.1 基本 40 → 40 控制。"}


def _duct_side_margin():
    """duct_layout 側向額外餘裕：參考橋 8 組 φ100、腹板 350、保護層 40。"""
    r = {m: duct_layout(8, 2, sec.yb - 1109, y_b=sec.yb, h=2100, side_margin=m) for m in (0, 15, 30)}
    return {"m0_cover_side": round(r[0].cover_side, 2), "m0_s_h": round(r[0].s_h, 2),
            "m15_cover_side": round(r[15].cover_side, 2), "m15_s_h": round(r[15].s_h, 2), "m15_ncol": r[15].n_col,
            "m30_ncol": r[30].n_col, "m30_e_max": round(r[30].e_max, 1), "m30_fits": r[30].fits,
            "_note": "預設排列把剩餘寬度全給列距，最外側管側向保護層恰 40(餘裕 0)；兩列配置只容 15 mm 額外餘裕，30 mm 退單列 e_max 1,029 < 1,109 排不下。"}


def _creep_aashto_A1():
    """AASHTO 5.4.2.3.2 潛變係數＋嚴格 AAEM λ（40 m 參考箱梁：V/S 箱內周長計 50%、H 75%、f'ci 32、t0 7 天）。"""
    import math
    vs = box_volume_surface(11000, 250, 5800, 200, 350, 2, 2100)
    c = aashto_creep(math.inf, 7, 75, vs, 32)
    rows = timing_sensitivity_aashto(7, [28, 90, 180], 0.0, -5.065 * 24.5 * 1600 / 8, 75, vs, 32)
    r28, r90 = rows[0], rows[1]
    return {"VS_mm": round(vs, 2), "psi_inf_t0_7": round(c.psi, 4),
            "t28_dphi": round(r28.dphi, 4), "t28_phi_r": round(r28.phi_r, 4),
            "t28_lam_exact": round(r28.lam_exact, 4), "t28_lam_approx": round(r28.lam_approx, 4),
            "t28_M_pier_kNm": round(r28.M_pier, 1),
            "t90_lam_exact": round(r90.lam_exact, 4), "t90_M_pier_kNm": round(r90.M_pier, 1),
            "t180_lam_exact": round(rows[2].lam_exact, 4),
            "_note": "台灣規範無 φ(t) 模式→採 AASHTO。嚴格 AAEM：分子 Δφ＝ψ(∞,t0)−ψ(t1−t0,t0)、分母 1+χψ(∞,t1)；"
                     "單一 Δφ 近似使 λ 偏大(28 天 0.513 vs 0.462)，對墩頂負彎矩與正束制偏保守。"
                     "ψ(∞,7)=1.30 遠小於 H7 示範的 φ∞=2.0。"}


def _duct_bundle_825_3():
    """§8.25.3 套管捆紮（每束 ≤3、端部 90 cm 內除外）：參考橋 8 組 φ100、保護層 40。"""
    b350, c350 = duct_layout_bundled(8, 2, sec.yb - 1109, web_t=350, y_b=sec.yb, h=2100)
    b300, c300 = duct_layout_bundled(8, 2, sec.yb - 1109, web_t=300, y_b=sec.yb, h=2100)
    return {"w350_none_e_max": round(c350["none"].e_max, 1), "w350_best": b350.bundle,
            "w350_best_e_max": round(b350.e_max, 1),
            "w300_none_e_max": round(c300["none"].e_max, 1), "w300_none_fits": c300["none"].fits,
            "w300_best": b300.bundle, "w300_best_e_max": round(b300.e_max, 1), "w300_best_fits": b300.fits,
            "_note": "腹板 350：垂直捆 e_max 1,169→1,189；腹板 300：不捆單列 e_max 1,029 排不下，水平捆兩列 1,169 排得下。"
                     "僅幾何可行性——腹板最小厚度、剪力、端部 90 cm 須展開之錨碇區未檢核。"}


def _duct_end_zone_825_3():
    """§8.25.3 端部 90 cm 內須回到 §8.25.2 間距：簡支拋物線 e_end=0、x≤900 取樣。"""
    f = parabolic_e(1109, 0, 40)
    es = [f(x / 1000) for x in range(0, 901, 100)]
    r350 = end_zone_duct_check(8, 2, es, sec.yb, 2100, web_t=350, bundle="V")
    r300 = end_zone_duct_check(8, 2, es, sec.yb, 2100, web_t=300, bundle="H")
    return {"e_hi_mm": round(r350.e_hi, 1),
            "w350_fits": r350.fits, "w350_layout": f"{r350.lay_hi.n_col}x{r350.lay_hi.n_row}",
            "w350_dx_max": r350.dx_max, "w350_dy_max": r350.dy_max,
            "w300_fits": r300.fits, "w300_layout": f"{r300.lay_hi.n_col}x{r300.lay_hi.n_row}",
            "w300_dx_max": r300.dx_max, "w300_dy_max": r300.dy_max,
            "_note": "端部 x≤0.9 m 偏心 ≤ 97.6，不捆紮排得下；過渡段腹板 300 水平捆→單列展開須橫移 50、豎移 140。"}


def _web_duct_A2():
    """腹板厚與管徑：台灣規範未規定；AASHTO 5.4.6.2（0.4）、C5.14.1.5.1c（最小腹板）、5.8.2.9（b_v 扣除）參考。"""
    a = web_duct_check(350, 100, 2, 2100)
    b = web_duct_check(300, 100, 2, 2100, bundle="H", dt_change=250)
    return {"w350_ratio": round(a.ratio, 4), "w350_ratio_ok": a.ratio_ok, "w350_web_min": a.web_min,
            "w350_bv_grouted": a.bv_grouted, "w350_bv_ungrouted": a.bv_ungrouted,
            "w300H_ratio": round(b.ratio, 4), "w300H_bundle_ratio": round(b.bundle_ratio, 4),
            "w300H_web_min_ok": b.web_min_ok, "w300H_bv_grouted": b.bv_grouted, "taper_250_mm": b.taper_min,
            "_note": "台灣三項均未規定（§8.9.3 僅漸變 12 倍）；以上為 AASHTO 參考值。H 捆紮束寬比 0.67 > 0.4，5.4.6.2 未定義捆紮 size，不判定。"}


def _deviation_force_A4():
    """過渡段偏折力（AASHTO 5.10.4.3）：腹板 300 水平捆→展開，橫移 50／豎移 140。"""
    Pu = 1.2 * 1395 * 19 * 140                       # 3.4.3：1.2 × 施拉力（19 股）
    r3 = deviation_force_check(Pu, 50, 140, 3000, 32, cover_side=40, cover_face=40)
    r6 = deviation_force_check(Pu, 50, 140, 6000, 32, cover_side=40, cover_face=40)
    return {"Pu_kN": round(Pu / 1e3, 1), "Lt3_R_lat_m": round(r3.R_lat / 1000, 1),
            "Lt3_F_out": round(r3.F_out, 1), "Lt3_F_in": round(r3.F_in, 1),
            "Vr_lat": round(r3.Vr_lat, 1), "Lt3_ok_lat": r3.ok_lat, "Lt3_ok_vert": r3.ok_vert,
            "Lt6_ok_lat": r6.ok_lat, "Lt6_ok_vert": r6.ok_vert,
            "Lt_min_lat_mm": round(r3.Lt_min_lat), "Lt_min_vert_mm": round(r3.Lt_min_vert),
            "s_max_mm": r3.s_max,
            "_note": "R=Lt²/(4Δ)（兩段反曲拋物線）。保護層 40 時，僅靠保護層抗拉脫需 Lt ≥ 4.06 m（垂直控制）；"
                     "較短者須配 tie-back／局部圍束。台灣規範未規定偏折力（§8.21.3 6.(2) 僅要求隔梁錨碇處配筋抵抗曲率分力）。"}


def _adjacent_duct_A7():
    """相鄰疊放曲線管道互推（AASHTO 5.10.4.3.1 三擇一）：跨中大半徑 vs 過渡段小半徑。"""
    Pu, Psv = 1.2 * 1395 * 19 * 140, 1115.0 * 19 * 140
    mid = adjacent_duct_radial_check(Pu, Psv, 180300, 2, 40, 32)
    trn = adjacent_duct_radial_check(Pu, Psv, 26500, 2, 40, 32)
    return {"F_each_mid": round(mid.F_each, 1), "F_stack_mid": round(mid.F_stack, 1),
            "Vr_between": round(mid.Vr_between, 1), "mid_spacing_ok": mid.spacing_ok,
            "F_each_trn": round(trn.F_each, 1), "trn_spacing_ok": trn.spacing_ok,
            "trn_s_req_mm": round(trn.s_req, 1), "trn_As_per_tie_mm2": round(trn.As_per_tie, 1),
            "tie_s_max": trn.tie_s_max,
            "_note": "跨中 R=180 m：F 24.7 ≤ V_r 151，淨距 40 已足；過渡段 R=26.5 m：F 168 > 151，"
                     "須淨距 ≥50 mm 或每 300 mm 配 267 mm² 圍束筋（f_s ≤ 0.6f_y），或內側管先灌漿再拉外側。"}


def _slab_thickness_TW():
    """台灣 §8.9.1／§8.9.2 箱梁翼板最小厚度；對照 RC 第七章 §7.1.22 底板 1/16。"""
    r = box_slab_thickness_TW(250, 200, 2400)
    w = box_slab_thickness_TW(250, 200, 6000)
    return {"top_req_2400": r.top_req, "bot_req_2400": r.bot_req, "rc_bot_req_2400": r.rc_bot_req,
            "top_req_6000": w.top_req, "bot_req_6000": w.bot_req, "bot_ok_6000": w.bot_ok,
            "_note": "PC 頂底板同為淨距/30（下限 15/14 cm，預鑄先拉 14/13）；RC 箱梁底板才是淨距/16。"
                     "淨距 6.0 m 時 PC 底板需 200 恰可，RC 規則會誤要求 375。"}


def _blister_stm_B1():
    """齒塊 STM 壓桿與節點（AASHTO 5.6.3 應變相容式 vs ACI β 表並列）：對齊算例_外置PT補強_錨固齒塊。"""
    r = blister_stm(3000, 3000, a_plate=200, b_plate=200, w_tie=150, theta_deg=45,
                    web_t=200, fc=40, node_type="CCC", A_node=400 * 200)
    r45 = blister_stm(3000, 3000, theta_deg=30, web_t=200, fc=40, node_type="CCC", A_node=400 * 200)
    return {"w_s_mm": round(r.w_s, 1), "A_cs_mm2": round(r.A_cs),
            "fcu_aashto_MPa": round(r.fcu_aashto, 2), "phiFns_aashto_kN": round(r.phiFns_aashto / 1e3),
            "fcu_aci_MPa": round(r.fcu_aci, 2), "phiFns_aci_kN": round(r.phiFns_aci / 1e3),
            "phiFnn_CCC_kN": round(r.phiFnn_aashto / 1e3), "node_ok": r.node_ok_aashto,
            "A_n_req_mm2": round(r.A_n_req), "A_cs_req_mm2": round(r.A_cs_req),
            "eps1": round(r.eps1, 5), "fcu_theta30_MPa": round(r45.fcu_aashto, 2),
            "node_CCT_aashto_kN": round(node_capacity_aashto(40, "CCT", 90000) / 1e3),
            "node_CCT_aci_kN": round(node_capacity(40, "CCT", 90000) / 1e3),
            "_note": "AASHTO 2008 無 β_s 表：f_cu=f'c/(0.8+170ε₁)≤0.85f'c、節點 0.85/0.75/0.65·φf'c；"
                     "ACI β 表(1.00/0.80/0.60)僅 CCC 與 AASHTO 同值，CCT 低 10%。"
                     "算例齒塊 400×250 於 200 腹板：壓桿與根部節點面積皆不足（與介面面積上限同一結論——須加大齒塊或改連續肋）。"}


def _pos_conn_detail_B2():
    """正彎矩接頭之錨定與配置（AASHTO 5.14.1.4.9b/c/d）：用 S4 之 M_req 設計兩種接頭。"""
    M = 5912.8                      # kN·m，0.6M_cr 控制（見 pos_moment_conn_S4）
    d = 1900.0
    st = pos_conn_strand(M, d, l_dsh=900, A_strand=140.0, projection=250)
    rb = pos_conn_rebar(M, d, ld=1200, dev_available=1500, bar_area=387.0)
    f900 = strand_stress_extended(900)
    return {"f_psl_600": round(strand_stress_extended(600)[0], 1),
            "f_pul_600": round(strand_stress_extended(600)[1], 1),
            "f_psl_900": round(f900[0], 1), "f_pul_900": round(f900[1], 1),
            "l_dsh_for_930MPa": round(strand_length_required(930), 1),
            "strand_n_req": round(st.n_req, 2), "strand_n_use": st.n_use,
            "strand_phiMn_kNm": round(st.phiMn), "strand_ok": st.ok,
            "rebar_n_req": round(rb.n_req, 2), "rebar_n_use": rb.n_use,
            "rebar_phiMn_kNm": round(rb.phiMn), "rebar_dev_ok": rb.dev_ok,
            "_note": "f_psl=(ℓ_dsh−203)/0.840、f_pul=(ℓ_dsh−203)/0.600（Eq.5.14.1.4.9c-1/-2，**式中無 d_b**）。"
                     "0.6M_cr 控制（5,913 kN·m）下：D22 需 24 支、φ15.2 延伸鋼絞線（ℓ_dsh 900）需 24 股——"
                     "兩者支數接近是巧合（鋼絞線應力高 2.8 倍、面積小 2.8 倍）。實際正束制僅 263 kN·m，最小量控制。"}


def _neg_conn_B5():
    """負彎矩接頭（5.14.1.4.8＋5.11.1.2.3＋5.14.1.4.6）：參考橋 40+40，墩頂負彎矩取包絡值。"""
    Mn = 13563.0            # kN·m，⑦ 兩等跨墩頂 Strength I 負彎矩（含體系轉換）
    r = negative_moment_connection(Mn, 2000, 40000, bar_area=387.0, db=22.2,
                                   embed_beyond_PI=3000, sigma_top=3.2, fc=40)
    nd = negative_moment_connection(Mn, 2000, 40000, embed_beyond_PI=1000, composite_deck=False)
    return {"As_req_mm2": round(r.As_req), "n_req": round(r.n_req, 2), "n_use": r.n_use,
            "one_third_bars": r.n_one_third, "embed_req_mm": r.embed_req,
            "embed_ok": r.embed_ok, "sigma_limit_bonded_MPa": round(r.sigma_limit, 2),
            "sigma_limit_unbonded_MPa": round(neg_top_tension_limit(40, False), 2),
            "service_ok": r.service_ok, "no_deck_embed_ok": nd.embed_ok,
            "no_deck_connection_required": nd.connection_required,
            "_note": "5.11.1.2.3：至少 1/3 負彎矩鋼筋延伸超過反曲點 ≥ max(d, 12d_b, 0.0625·淨跨)；"
                     "本例 0.0625×40,000＝2,500 控制（d 2,000、12d_b 266）。"
                     "5.14.1.4.6 梁頂拉應力用 Table 5.9.4.1.2-1 但以 f'c 代 f'ci → 有握裹 0.63√f'c＝3.98 MPa。"}


def _multi_stage_C1():
    """多階段逐跨施工：三跨 40 m，每跨 120 天，各跨材齡不同 → 各階段自有 λ。"""
    spans, wdl = [40.0, 40.0, 40.0], 40.4
    x = 40.0                                   # 第一中墩
    pairs = stage_moments_span_by_span(spans, wdl, x)
    sched = span_by_span_schedule(3, 120)
    st = [StageSpec(f"span{k+1}", t0, tc, mi, mii)
          for k, ((mi, mii), (t0, tc)) in enumerate(zip(pairs, sched))]
    r = multi_stage_redistribution(st)
    slow = multi_stage_redistribution(
        [StageSpec(f"span{k+1}", t0, tc, mi, mii)
         for k, ((mi, mii), (t0, tc)) in enumerate(zip(pairs, span_by_span_schedule(3, 400)))])
    return {"M_I_total": round(r.M_I_total, 1), "M_II_total": round(r.M_II_total, 1),
            "lam_per_stage": [round(row[1], 4) for row in r.rows],
            "M_per_stage": [round(row[4], 1) for row in r.rows],
            "M_total_multi": round(r.M_total, 1), "lam_equiv": round(r.lam_equiv, 4),
            "M_single_lam": round(r.M_single, 1), "single_lam": round(r.single_lam, 4),
            "age_gap_days": r.age_gap, "age_gap_warn": r.age_gap_warn,
            "slow_age_gap": slow.age_gap, "slow_warn": slow.age_gap_warn,
            "slow_M_total": round(slow.M_total, 1),
            "_note": "逐跨施工各跨材齡不同→λ 各異（先建的跨剩餘潛變少、λ 小：0.099／0.167／0.586）。"
                     "🔴 本例單一 λ 給 −640、逐階段疊加只有 −335——但反向案例（各階段同號）單一 λ 反而偏小，"
                     "**方向依各階段 M_I／M_II 的正負組合而定，不可假設簡化偏保守**（見測試 test_multi_stage_C1）。"
                     "每跨 400 天時材齡差 800 天 > 365，便覽 §3.5.4 但書要求分步考慮施工中斷面力與材齡差。"}

# ── C2：懸臂工法 X₁／X₀ 自動組成（H3 幾何 → H7 重分配）────────────────────
# 參考橋＝算例_懸臂工法（80+80 m，中墩出平衡懸臂 ±40 m，兩端 40 m 支架段）。
_H3_SEG = [("0號塊", 2.0, 605.0), ("S1", 6.75, 643.0), ("S2", 11.25, 599.0),
           ("S3", 15.75, 550.0), ("S4", 20.25, 497.0), ("S5", 24.75, 439.0),
           ("S6", 29.25, 385.0), ("S7", 33.75, 353.0), ("S8", 38.25, 341.0)]


def _cantilever_C2():
    spans, piers, closures = [80.0, 80.0], [80.0], [40.0, 120.0]
    units = cantilever_units(spans, piers, closures)
    arm_kN = sum(g for _, _, g in _H3_SEG)                 # 單臂 4,412 kN
    w_eq = 2 * arm_kN / 80.0                               # 兩臂等效 110.3 kN/m
    r = cantilever_X1X0(spans, units, w_eq, dphi=1.0, w_sdl=15.0)
    ip = min(range(len(r.xs)), key=lambda k: abs(r.xs[k] - 80.0))
    # H3 交叉驗證：逐節塊集中載重 + 掛籃（尖端外 0.5 m）→ 0 號塊端 x=4 m 斷面
    lo = [(80.0 + d, g) for _, d, g in _H3_SEG] + [(120.5, 800.0)]
    M_cant = cantilever_M_I_points(units, lo, 84.0, "R", overhang=1.0)
    # 逐節塊多階段（各節塊材齡不同）：10 天/節塊、合龍再 30 天
    # 兩臂同步推進（各節塊左右同齡）；M_I 於墩頂取右側自由體＝只有右臂載重，
    # M_II 則兩臂皆計——同一組載重的兩個平衡解，故 ΔM 仍為支承間線性。
    segs, sch = [], []
    base = cantilever_schedule(9, 10.0, t0=7.0, days_to_closure=30.0)
    for (nm, d, g), ts in zip(_H3_SEG, base):
        segs += [(nm + "R", 80.0 + d, g), (nm + "L", 80.0 - d, g)]
        sch += [ts, ts]
    st = cantilever_stages(spans, units, segs, 80.0, sch)
    ms = multi_stage_redistribution(st)
    dM = lambda xx: sum(sp.M_II - sp.M_I for sp in
                        cantilever_stages(spans, units, segs, xx, sch))
    d0, d80 = dM(0.0), dM(80.0)
    lin_err = max(abs(dM(xx) - (d0 + (d80 - d0) * xx / 80.0))
                  for xx in (10.0, 20.0, 40.0, 60.0, 70.0))
    return {
        "config": "80+80m 中墩平衡懸臂±40m/合龍40,120/支架段2×40m/w_eq110.3/SDL15/Δφ=1.0",
        "arms": [[u.a_left, u.a_right] for u in units],
        "layout_3span": [list(cantilever_layout([40.0, 50.0, 40.0])[0]),
                         list(cantilever_layout([40.0, 50.0, 40.0])[1])],
        "layout_4span_arms": [[u.a_left, u.a_right] for u in cantilever_units(
            [60.0, 100.0, 100.0, 60.0], *cantilever_layout([60.0, 100.0, 100.0, 60.0]))],
        "falsework": [[round(a, 1), round(b, 1)] for a, b in cantilever_falsework(spans, units)],
        "w_eq_kNpm": round(w_eq, 2),
        "lam": round(r.lam, 4),
        "X1_pier": round(r.piers[0][1], 1), "X0_pier": round(r.piers[0][2], 1),
        "Xfinal_pier": round(r.piers[0][3], 1),
        "X1_closure": round(r.closures[0][1], 1) + 0.0, "X0_closure": round(r.closures[0][2], 1),
        "Xfinal_closure": round(r.closures[0][3], 1),
        "unbalanced_pier": round(r.unbalanced[0][1], 1) + 0.0,
        "linear_ok": r.linear_ok,
        "M_post_pier": round(r.M_post[ip], 1),
        "M_dead_total_pier": round(r.M_dead_total[ip], 1),
        "M_lumped_wrong_pier": round(r.M_lumped_wrong[ip], 1),
        "lumped_err_pier": round(r.lumped_err_pier, 1),
        "M_cant_x4_with_traveler": round(M_cant, 2),
        "stage_lam": [round(row[1], 4) for row in ms.rows[::2]],
        "stage_dM_linear_err": round(lin_err, 9),
        "stage_M_I_total": round(ms.M_I_total, 1),
        "stage_M_II_total": round(ms.M_II_total, 1),
        "stage_M_total": round(ms.M_total, 1),
        "stage_lam_equiv": round(ms.lam_equiv, 4),
        "_note": "X₁＝懸臂體系（合龍前瞬間）、X₀＝同載重一次完工；X_final=X₁+λ(X₀−X₁)。"
                 "墩頂 X₁=−88,240＝−w·a²/2 閉合解；合龍點 X₁=0（自由端）→ 潛變生成正束制彎矩。"
                 "🔴 M_cant_x4_with_traveler=−90,861.25 與 construction.cantilever_moment（手填力臂）"
                 "逐位元相同——H3 幾何自動組成與既有 H3 算路互為交叉驗證。"
                 "🔴 lumped_err_pier：把合龍後載重（支架段落架＋SDL）一併塞進 X₀ 再乘 λ，"
                 "墩頂負彎矩少算 22,491 kN·m（19.2%）、**偏不安全**——X₀ 與 X₁ 必須是同一組載重。"
                 "stage_dM_linear_err≈0 為逐節塊版之自我驗證（束制場自平衡→支承間線性）。"}


golden = {
    "_about": "40m參考橋黃金答案(台灣HS20-44/2車道/8組×19股最小設計)。Python引擎與JS網頁前端共用驗證源。由 make_golden.py 自動產生，請勿手改。",
    "influence_simple_40m": {
        "peak_M_mid_m": round(il_moment_peak(40, 20), 2),
        "peak_M_a10_m": round(il_moment_peak(40, 10), 2),
        "truck_absmax_kNm": round(taiwan_truck_moment(40)),
        "lane_M_kNm": round(taiwan_lane_moment(40)),
        "per_lane_M_LL_IM_kNm": round(taiwan_per_lane_moment(40)),
        "envelope_peak_kNm": round(max(m for _, m in taiwan_moment_envelope(40))),
        "envelope_at_L4_kNm": round(dict((round(a), m) for a, m in taiwan_moment_envelope(40))[10]),
    },
    "loads": {
        "M_LL_IM_2lane_kNm": round(M_LL),
        "StrengthI_kNm": round(c["Strength_I"]),
        "ServiceI_kNm": round(c["Service_I"]),
        "ServiceIII_kNm": round(c["Service_III"]),
    },
    "prestress": {
        "config": "8組×19股", "Aps_mm2": round(ten.Aps),
        "loss_pct": round(L.loss_pct * 100, 1), "fpe_MPa": round(L.fpe),
        "Pe_kN": round(L.Pe / 1e3),
        "Pe_min_kN": round(Pe_min_zero_tension(sec, ten.e, c["Service_I"]) / 1e3),
    },
    "service": {"sigma_bot_MPa": round(sb, 2), "sigma_top_MPa": round(st, 2)},
    "shear_D1": {"fpc_MPa": round(sh.fpc, 2), "Vp_kN": round(sh.Vp / 1e3),
                 "sigma1_MPa": round(sh.sigma1, 2), "Vcw_kN": round(sh.Vcw / 1e3),
                 "Vs_req_kN": round(sh.Vs_req / 1e3)},
    "flexure_M1": {"c_mm": round(fx.c, 1), "fps_MPa": round(fx.fps),
                   "Mn_kNm": round(fx.Mn), "CR": round(fx.CR, 2)},
    "deflection": {"LBR_pct": round(df.LBR * 100, 1), "delta_LL_mm": round(df.d_LL, 1),
                   "net_longterm_mm": round(df.net_long_term, 1), "camber_mm": round(df.camber)},
    "fatigue_P1": (lambda fa: {"dsig_ps_MPa": round(fa.dsig_ps, 1),
                               "sig_c_max_MPa": round(fa.sig_c_max, 2),
                               "stirrup_250_MPa": round(stirrup_fatigue(565, 250, 402, 1692)[0]),
                               "stirrup_150_MPa": round(stirrup_fatigue(565, 150, 402, 1692)[0])})
                  (fatigue_check(sec, L.Pe, ten.e, 28800, 3222, 40)),
    "torsion_D2": (lambda tr: {"fpc_MPa": round(tr.fpc, 2), "Tcr_kNm": round(tr.Tcr),
                               "threshold_kNm": round(tr.threshold), "neglect_explicit": tr.neglect})
                  (torsion_check(sec, L.Pe, 40, 23.1e6, 26200, 1900)),
    "transverse_D3": {"cantilever_phiMn_kNm": round(slab_flexure(105.8, 1571, 200, 40, 420).phiMn, 1),
                      "span_phiMn_kNm": round(slab_flexure(133.8, 2172, 200, 40, 420).phiMn, 1),
                      "support_phiMn_kNm": round(slab_flexure(150.3, 2534, 200, 40, 420).phiMn, 1),
                      "As_min_mm2": round(As_min_slab(40, 420, 1000, 200))},
    "temperature_T1": temp_gradient_AASHTO(18.0, 5.0, True),
    "bearing_E1": (lambda b: {"R_max_kN": round(1440+R_LL_HS20), "gamma_s": round(b.gamma_s,2),
                              "shape_S": round(b.shape_S,1), "sigma_TL_MPa": round(b.sigma_TL,2),
                              "sigma_TL_limit_MPa": round(b.sigma_TL_limit,2), "H_m_kN": round(b.H_m,1),
                              "gamma_ok": b.gamma_ok, "sigma_ok": b.sigma_ok,
                              "stability_ok": b.stability_ok, "H_ok": b.H_ok, "no_uplift": b.no_uplift})
                  (bearing_check(1440+R_LL_HS20, 1440, R_LL_HS20, 40, 100, 550, 450, te=10, G_kgf=8)),
    "anchorage_F1": {"Pu_kN": round(an.Pu), "sum_Tburst_kN": round(an.sum_Tburst),
                     "Fspall_kN": round(an.Fspall), "As_spall_mm2": round(an.As_spall),
                     "spiral_Pult_kN": round(sp[0]), "bearing_margin": round(sp[1], 2),
                     "bearing_ok": sp[2]},
    "expansion_E2": (lambda j: {"shortening_mm": round(j.shortening,1), "g_max_mm": round(j.g_max,1),
                                "capacity_mm": j.capacity, "joint": j.joint_type})(expansion_joint(8.8,12.6,8.0,20)),
    "live_load_TW_HS20_40m": {"model": "卡車或車道取大", "impact": round(taiwan_impact(40),4),
        "per_lane_M_kNm": round(taiwan_per_lane_moment(40)), "per_lane_V_kN": round(taiwan_per_lane_shear(40)),
        "truck_V_support_kN": round(taiwan_truck_shear(40), 1), "per_lane_V_30m_kN": round(taiwan_per_lane_shear(30), 1),
        "_note": "台灣 HS20-44；對照 HL-93 每車道 5673/588。卡車雙向掃描：40m 卡車 301.1 < 車道 304（車道控制）；30m 卡車控制"},
    # 標準算例 T1（自含斷面 h=2000、無預力 σ_base=1.2）：孤島保守 illustration → +2.00 控制
    "temperature_SE_T1": (lambda r: {"Tu_C": round(r.Tu,2), "TL_C": round(r.TL,1),
        "sigSE_top_MPa": round(r.sigma_pos["頂板頂"],2), "sigSE_bot_pos_MPa": round(r.sigma_pos["底板底"],2),
        "sigSE_bot_neg_MPa": round(r.sigma_neg["底板底"],2),
        "service_bot_neg_MPa": round(thermal_service_check(r.sigma_neg["底板底"],1.2,0.5)[0],2),
        "service_ok": thermal_service_check(r.sigma_neg["底板底"],1.2,0.5)[1],
        "_note": "自含斷面 h=2000、σ_base=1.2（無預力）孤島 illustration；實際橋見 temperature_integrated_T1"})(
        self_equilibrating_stress(
            [ThermalBand(0,300,3_000_000,11.5),ThermalBand(300,400,80_000,2.5),
             ThermalBand(400,1750,1_080_000,0),ThermalBand(1750,2000,1_375_000,0)],
            1.26e12,870,2000,[("頂板頂",0,18.0),("底板底",2000,0.0)])),
    # ★ 接線：config A 斷面 + 引擎實際服務性底緣（含預力）→ 真實參考橋的 T1 整合檢核
    "continuous_pier": _continuous_pier(),
    "cont_tendon_force_default": _cont_tendon_force_default(),
    "cont_envelope_taiwan": _cont_envelope_taiwan(),
    "pier_cap_tendon": _pier_cap_tendon(),
    "variable_section_haunch": _variable_section_haunch(),
    "cont_shear_taiwan": _cont_shear_taiwan(),
    "rear_axle_scan": _rear_axle_scan(),
    "simple_shear_dv": _simple_shear_dv(),
    "pier_cap_tendon_losses": _pier_cap_tendon_losses(),
    "cont_shear_scan": _cont_shear_scan(),
    "tendon_slip_G1": _tendon_slip_G1(),
    "mid_anchor_G1": _mid_anchor_G1(),
    "loss_profile_G1": _loss_profile_G1(),
    "blister_B1": _blister_B1(),
    "staging_redist_S1": _staging_redist_S1(),
    "staged_envelope_S2": _staged_envelope_S2(),
    "duct_size_T83": _duct_size_T83(),
    "cont_envelope_gamma_min": _cont_envelope_gamma_min(),
    "staged_shear_S3": _staged_shear_S3(),
    "pos_conn_detail_B2": _pos_conn_detail_B2(),
    "neg_conn_B5": _neg_conn_B5(),
    "multi_stage_C1": _multi_stage_C1(),
    "pos_moment_conn_S4": _pos_moment_conn_S4(),
    "durability_cover_T12": _durability_cover_T12(),
    "duct_side_margin": _duct_side_margin(),
    "creep_aashto_A1": _creep_aashto_A1(),
    "duct_bundle_825_3": _duct_bundle_825_3(),
    "duct_end_zone_825_3": _duct_end_zone_825_3(),
    "web_duct_A2": _web_duct_A2(),
    "deviation_force_A4": _deviation_force_A4(),
    "adjacent_duct_A7": _adjacent_duct_A7(),
    "blister_stm_B1": _blister_stm_B1(),
    "slab_thickness_TW": _slab_thickness_TW(),
    "temperature_integrated_T1": (lambda r: {"section": "配置A h=2100", "Tu_C": round(r.Tu,2), "TL_C": round(r.TL,2),
        "sigSE_bot_neg_MPa": round(r.sigma_neg["底板底"],2), "service_base_MPa": round(sb,2),
        "service_total_MPa": round(thermal_service_check(r.sigma_neg["底板底"], sb, 0.5)[0],2),
        "service_ok": thermal_service_check(r.sigma_neg["底板底"], sb, 0.5)[1],
        "_note": "σ_base=引擎服務性底緣(含預力,HS20)；預力餘裕吸收熱應力→通過(孤島illustration的+2.00為自含斷面/無預力假象)"})(
        self_equilibrating_stress(
            [ThermalBand(0,250,11000*250,12.58),ThermalBand(250,300,700*50,6.08),
             ThermalBand(300,400,700*100,2.5),ThermalBand(400,1900,700*1500,0),
             ThermalBand(1900,2100,5800*200,0)],
            sec.I, sec.h-sec.yb, sec.h, [("底板底",2100,0.0)], Ec=29700)),
    "tendon_profile_G1": {
        "config": "8組×21股 (HL-93 KB敘述軌)", "a_mm": 1109,
        "theta_end_rad": round(g1.theta_end, 3), "R_m": round(g1.R / 1000, 1),
        "w_eq_transfer_kNpm": round(g1.w_eq_transfer, 1), "w_eq_service_kNpm": round(g1.w_eq_service, 1),
        "w_DL_kNpm": round(w_DL_g1, 1), "LBR_transfer": round(g1.LBR_transfer, 3),
        "LBR_service": round(g1.LBR_service, 3), "friction_single_end": round(g1.fric_single_end, 3),
        "friction_dual_mid": round(g1.fric_dual_mid, 3), "R_ok": g1.R_ok,
        "_note": "拋物線等效荷載法 w_eq=8Pa/L²；對齊算例_鋼腱線形設計(8組×21股HL-93軌)；w_DL=8(M_DC+M_DW)/L²=144與參考橋自洽"},
    "duct_layout_G1": {
        "config": "參考橋 8組/2腹板/φ100管/腹板350/保護層40/垂直淨距40/骨材25/rule=tw(台灣明文)",
        "n_col": dl.n_col, "n_row": dl.n_row,
        "s_req_mm": round(dl.s_req, 1), "s_h_mm": round(dl.s_h, 1),
        "web_avail_mm": round(dl.web_avail, 1), "pitch_v_mm": round(dl.pitch_v, 1),
        "rows_y_mm": [round(y, 1) for y, _ in dl.rows],
        "y_cgs_mm": round(dl.y_cgs, 1), "y_bot_mm": round(dl.y_bot, 1),
        "cover_bot_mm": round(dl.cover_bot, 1),
        "cover_ok": dl.cover_ok, "s_v_ok": dl.s_v_ok, "s_h_ok": dl.s_h_ok, "fits": dl.fits,
        "e_max_mm": round(dl.e_max, 1), "e_max_point_mm": round(dl.e_max_point, 1),
        "_note": "實配排列：每腹板由下往上逐層填滿、整體平移使形心=分析用CGS。rule=tw 淨距需求 max(40,1.5d_agg)=40 → 兩列兩層，e=1109 排得下。e_max(實配) < e_max_point(單點簡化)"},
    "duct_layout_G1_strict": {
        "config": "同上但保守側：保護層75(PTBG實務)/垂直淨距100/rule=od(淨距再取≥孔道外徑)",
        "n_col": dls.n_col, "n_row": dls.n_row,
        "s_req_mm": round(dls.s_req, 1), "web_avail_mm": round(dls.web_avail, 1),
        "pitch_v_mm": round(dls.pitch_v, 1), "rows_y_mm": [round(y, 1) for y, _ in dls.rows],
        "cover_bot_mm": round(dls.cover_bot, 1), "cover_ok": dls.cover_ok, "fits": dls.fits,
        "e_max_mm": round(dls.e_max, 1),
        "_note": "同一配置在保守側排不下(單列四層、最底層管外緣落在梁底以下130)。兩案並列以顯示淨距規則的影響：φ100時 rule=od 使需求由40變100，直接決定腹板內排一列或兩列"},
    "duct_layout_G1_pier": {
        "config": "墩頂斷面 e_pier=-900(形心上方)/保護層40/垂直淨距40/rule=tw",
        "n_col": dlp.n_col, "n_row": dlp.n_row,
        "y_cgs_mm": round(dlp.y_cgs, 1), "y_top_mm": round(dlp.y_top, 1),
        "top_ok": dlp.top_ok, "fits": dlp.fits,
        "e_pier_limit_mm": round(sec.yb - 2100 + 40 + 50 + (dlp.y_top - dlp.y_cgs), 1),
        "_note": "墩頂腱在形心上方，控制的是頂緣保護層。e_pier=-900 使頂層管外緣穿出梁頂；可行下限 = yb - h + cover + od/2 + (y_top - y_cgs)"},
    "friction_profile_G1": {
        "config": "參考橋 L=40m / a=1109 / mu=0.25 / K=0.003",
        "kappa_rad_per_m": round(fsg[0][2], 6),
        "alpha_mid_rad": round(friction_angle(fsg, 20, 40), 5),
        "alpha_full_rad": round(friction_angle(fsg, 40, 40), 5),
        "both":  {"start": round(fps["both"].at_start, 4),  "mid": round(fps["both"].at_mid, 4),
                  "end": round(fps["both"].at_end, 4),      "avg": round(fps["both"].avg, 4),
                  "max": round(fps["both"].max_ratio, 4),   "x_max_m": round(fps["both"].x_max, 1)},
        "start": {"start": round(fps["start"].at_start, 4), "mid": round(fps["start"].at_mid, 4),
                  "end": round(fps["start"].at_end, 4),     "avg": round(fps["start"].avg, 4),
                  "max": round(fps["start"].max_ratio, 4),  "x_max_m": round(fps["start"].x_max, 1)},
        "end":   {"start": round(fps["end"].at_start, 4),   "mid": round(fps["end"].at_mid, 4),
                  "end": round(fps["end"].at_end, 4),       "avg": round(fps["end"].avg, 4),
                  "max": round(fps["end"].max_ratio, 4),    "x_max_m": round(fps["end"].x_max, 1)},
        "alt":   {"start": round(fps["alt"].at_start, 4),   "mid": round(fps["alt"].at_mid, 4),
                  "end": round(fps["alt"].at_end, 4),       "avg": round(fps["alt"].avg, 4),
                  "max": round(fps["alt"].max_ratio, 4),    "x_max_m": round(fps["alt"].x_max, 1)},
        "_note": "1−e^−(μα+Kx)，α 自張拉端累積(∫|e″|dx，反曲須逐段取絕對值)。跨中四配置同值(左右對稱、路徑各半)=0.0840 與 tendon_profile_G1.friction_dual_mid 一致；單端遠端 0.1609 與 friction_single_end 一致。alt 為該斷面半數腱自左半數自右的平均，個別腱仍是單端分佈"},
    "tendon_forces_G1": {
        "config": "x=L/4=10m / fpj=1395 / Ap_each=140x19 / 其他損失 180 MPa / 交錯端",
        "single_col_4row": {"Pe_kN": round(tfA.Pe_total / 1e3, 1), "e_eff_mm": round(tfA.e_eff, 2),
                            "e_geom_mm": round(tfA.e_geom, 2), "de_mm": round(tfA.de, 2),
                            "x_eff_mm": round(tfA.x_eff, 2)},
        "two_col_2row":    {"Pe_kN": round(tfB.Pe_total / 1e3, 1), "e_eff_mm": round(tfB.e_eff, 2),
                            "e_geom_mm": round(tfB.e_geom, 2), "de_mm": round(tfB.de, 2),
                            "x_eff_mm": round(tfB.x_eff, 2)},
        "single_col_start":{"Pe_kN": round(tfS.Pe_total / 1e3, 1), "de_mm": round(tfS.de, 2),
                            "x_eff_mm": round(tfS.x_eff, 2)},
        "Pe_equals_avg_ratio": abs(tfA.Pe_total - tfA.Pe_avg_ratio) < 1.0,
        "_note": "總 Pe 以平均損失率算等價(各項對 f_pe 線性)，差別在合力位置。單列四層交錯時序號＝層號→低層自起點損失小→合力下移 Δe=+5.10；兩列兩層同層左右配對→垂直相消 Δe=0 但橫向 x_eff=-4.30；全部自起點則各腱損失相同→Δe=x_eff=0"},
    "stm_F2": {
        "config": "8組×19股 (參考橋, 端橫隔版 General Zone)",
        "sigma_pe_MPa": round(f2.sigma_pe, 2), "T_burst_kN": round(f2.T_burst / 1e3),
        "d_burst_mm": round(f2.d_burst), "As_burst_mm2": round(f2.As_burst),
        "beta_strut_required": f2.beta_strut_required, "fcu_node_CCT_MPa": round(f_cu(40, 0.80), 1),
        "node_A_CCT_phiFnn_kN": round(node_capacity(40, "CCT", 90000) / 1e3),
        "_note": "AASHTO 5.9.5.6.3 General Zone T_burst=0.25ΣP(1−a/h)；對齊算例_STM端橫隔版設計(參考橋8組×19股)"},
    "durability_N1": {
        "config": "standalone 卡；codified 驗收門檻/限值",
        "grout_w_c_max": GROUT["w_c_max"], "grout_f28_min_MPa": GROUT["f28_min"],
        "grout_bleed_max_pct": GROUT["bleed_max_pct"], "grout_chloride_max_pct": GROUT["chloride_max_pct"],
        "rebar_permanent_max_MPa": rebar_stress_limit("常時"),
        "rebar_fatigue_general_max_MPa": rebar_stress_limit("疲勞_一般"),
        "rebar_fatigue_deck_max_MPa": rebar_stress_limit("疲勞_床版翼緣"),
        "pc_fatigue_rule": "min(0.60Pu, 0.75Py)", "design_life_TW": list(design_life("台灣")),
        "design_life_AASHTO": design_life("AASHTO"), "design_life_JP": design_life("日本"),
        "_note": "灌漿驗收 道示Ⅲ17.6.6(2)；100年應力限制 道示Ⅲ6.2.2/6.3.2；混凝土疲勞壓應力(表-6.3.5)canonical在C1(2026-07-01 NLM查證:箱形60→18.0/80→26.0,N1原誤已校正)"},
    "construction_stage_H1H2": {
        "config": "40m參考橋 8組×19股，支架/托架上施拉 fci=32",
        "transfer_tension_limit_MPa": round(transfer_tension_limit(32), 2),
        "S2_fullPT_top_MPa": round(h8.sigma_top, 2), "S2_fullPT_top_ok": h8.top_ok,
        "S2_fullPT_bot_MPa": round(h8.sigma_bot, 2),
        "S2_batch4_top_MPa": round(h4.sigma_top, 2), "S2_batch4_top_ok": h4.top_ok,
        "S2_batch4_bot_MPa": round(h4.sigma_bot, 2),
        "S3_strike_top_MPa": round(h3s.sigma_top, 2), "S3_strike_bot_MPa": round(h3s.sigma_bot, 2),
        "_note": "H1/H2 支架上施拉自重未活化(M_sw=0)→過平衡頂緣引張；全PT頂拉+底壓(0.55f'ci=17.6,一般橋梁裁示)皆超限→分批4組通過(頂0.93/底-9.59)；脫架自重活化回壓。對齊算例_40m參考橋施工階段應力歷程"},
    "cantilever_H3": {
        "config": "80+80m 變深連續梁 h_pier4.5/h_mid2.2，8節塊/側",
        "h_mid_m": round(variable_depth(0, 4.5, 2.2, 40), 2), "h_pier_m": round(variable_depth(40, 4.5, 2.2, 40), 2),
        "h_at_x20_m": round(variable_depth(20, 4.5, 2.2, 40), 3),
        "M_selfweight_about_x4_kNm": round(cantilever_moment(w_h3, a_h3)),
        "M_cant_max_kNm": round(cantilever_moment(w_h3, a_h3, 800, 36.5)), "M_cant_max_published_kNm": 90825,
        "delta_elastic_mm": 147, "delta_long_term_mm": round(long_term_deflection(147, 2.0)),
        "_note": "變深h(x)端點2.2/4.5；懸臂彎矩Σ(G·arm)+掛籃；力臂統一對 x=4 斷面（掛籃距墩CL 40.5m→臂36.5m；2026-09-15 修正原混用參考點的94,025）；長期下撓δ(1+φ)。對齊算例_懸臂工法施工階段設計"},
    "launching_H4": {
        "config": "40m等跨 ILM 推進，等深 h=2.2m A=4.870m² Zb=3.093e9mm³",
        "M_cantilever_neg_kNm": round(launching_cantilever_moment(120.0, 14.0)),
        "M_span_pos_kNm": round(Mpos_h4), "Pc_required_kN": round(Pc_h4), "Pc_published_kN": 45100,
        "n_tendons": n_tendons(Pc_h4, 2510),
        "sigma_bot_min_MPa": round(launching_bottom_stress(45100, 4.870e6, Mpos_h4, 3.093e9), 2),
        "jacking_force_kN": round(jacking_force(0.10, 120 * 200)),
        "bearing_stress_widened_MPa": round(bearing_stress(7200, 960000), 2),
        "_note": "ILM 推進包絡 M⁻/M⁺；臨時置中預力(e=0)Pc=(M⁺/Zb+σ_res)·A→底緣恰餘1.5壓；頂推F=μ·W；支壓R/A。對齊算例_推進工法施工設計"},
    "segmental_H5H6": {
        "config": "H5 預鑄SBS 40m Ac=4.20m² + H6 接縫BCM 80m Ac=2.85m²",
        "H5_segment_weight_kN": round(segment_weight(4.20, 2.5, 25), 1),
        "H5_joint_min_prestress_kN": round(joint_min_prestress(4.20e6), 1),
        "H5_joint_compression_4tendon_MPa": round(joint_compression(4 * 480, 4.20e6), 3),
        "H6_joint_min_compression_kN": round(joint_min_prestress(2.85e6), 1),
        "H6_shear_key_design_capacity_kN": round(Vkey_h6, 1),
        "H6_LS3_utilization": round(shear_key_utilization(2850, 20, Vkey_h6), 2),
        "H6_bonded_pt_ratio": round(bonded_pt_ratio(14000, 42500), 3),
        "_note": "H5 節塊重Ac·L·γ/拼裝期接縫壓0.21；H6(道示法)剪力鍵V_fuk·ξ/LS3驗核比/黏結PT≥30%。對齊算例_預鑄節塊工法施工設計+算例_節塊接縫設計"},
    "design_inverse": {
        "config": "40m參考橋反解設計庫閉環（設計=驗算之逆）",
        "Pe_min_zero_tension_kN": round(Pe_min_ref / 1e3),
        "min_tendon_groups": d_ngroups,
        "required_drape_LBR0985_mm": round(d_drape),
        "min_Sb_zero_tension_e9mm3": round(d_Sbmin / 1e9, 3),
        "actual_Sb_e9mm3": round(sec.Sb / 1e9, 3),
        "_note": "階段4反解設計庫(純Python零相依,SymPy為推導輔助)：min_tendon_groups(Pe_min)=8=實際8組；required_drape(LBR0.985)=1109=e_m；min_Sb=1.992e9≤實際2.473e9(斷面足夠)。Pe_min與min_Sb為σ_b=0同式對偶。pint單位QA/SymPy驗證工具未安裝,屬開發期另檔,不進零相依runtime"},
    "seismic_S1": {
        "config": "落橋防止(台灣耐震§8.5) 40m跨/基面起H=10m/第二類地盤 Le=50m",
        "min_NL_cm": round(seis.min_falloff_length(40, 10, 0), 1),
        "min_NL_skew30_cm": round(seis.min_falloff_length(40, 10, 30), 2),
        "u_G_cm": round(seis.ground_relative_displacement("第二類", 5000, 1.2), 2),
        "L_N_required_cm": round(seis.required_falloff_length(70, 15, 22.5), 1),
        "restrainer_Fy_kN": round(seis.restrainer_yield_strength(1800), 1),
        "_note": "min N_L=(50+0.25L+H)(1+S²/8000)=70cm；u_G=ε_G·Le·(S_III/S_II)=0.00375·5000·1.2=22.5cm；活動支承L_N=max(70,15+22.5)=70cm(min N_L控制)；裝置F_y=1.5R_d。對齊算例_落橋防止系統設計"},
    "seismic_S2": {
        "config": "隔震LRB等效線性化迭代(台灣耐震第7章) W=8000kN/Q_d=400kN/K_d=6kN/mm/剛性墩/S_II,1=0.60",
        "D_d_mm": round(seis_S2.D_d * 1000, 1),
        "T_e_s": round(seis_S2.T_e, 3),
        "K_eff_kN_per_m": round(seis_S2.K_eff, 1),
        "xi_e_pct": round(seis_S2.xi_e * 100, 2),
        "B1": round(seis_S2.B1, 4),
        "V_b_kN": round(seis_S2.V_b_secant, 1),
        "iterations": seis_S2.iterations,
        "B1_15pct_table3_1": round(seis.damping_correction_B1(0.15), 4),
        "_note": "等效線性化迭代:s_D→D_d→K_eff/ξ_eq→T_e→B_1(表3-1內插)→S_a→新s_D。收斂D_d=221mm/T_e=2.03s/ξ_e=14.75%/V_b=1727kN(6次)。B_1(15%)=1.375為表3-1線性內插(算例舊用1.4近似)。對齊算例_隔震與消能設計"},
    "seismic_S3": {
        "config": "橋墩韌性容量設計(台灣耐震§4.2/5.3) 圓D150/fc280/fyh2800/Ag17671/Ac15394/Pe800000kgf/Mn3000tfm",
        "M_p_tfm": round(seis.overstrength_moment(3000), 1),
        "V_u_tf": round(seis.capacity_shear(3900, 8), 1),
        "rho_s_pct": round(seis.rho_s_circular(280, 2800, 17671, 15394, 800000) * 100, 3),
        "spacing_limit_cm": round(seis.confinement_spacing_limit(150, 3.6), 1),
        "plastic_hinge_L_cm": round(seis.plastic_hinge_length(150, 800), 1),
        "_note": "M_p=1.3Mn=3900；V_u=ΣM_p/L_c=487.5；圓柱ρ_s=max(式5-5=0.67%,式5-6=0.84%)=0.84%(軸力式控制)；間距≤min(15,D/4,6db)=15；ℓ0=max(柱深150,ℓc/6=133,45)=150。對齊算例_橋墩韌性耐震設計"},
    "seismic_S5": {
        "config": "液狀化土壤參數折減D_E(台灣耐震表8-1) 三級距×深度×R_s",
        "DE_L1_shallow_loose": round(seis.liquefaction_reduction_DE(0.3, 5, 0.2), 4),
        "DE_L1_shallow_dense": round(seis.liquefaction_reduction_DE(0.3, 5, 0.4), 4),
        "DE_L1_deep": round(seis.liquefaction_reduction_DE(0.3, 15, 0.2), 4),
        "DE_L2_shallow_loose": round(seis.liquefaction_reduction_DE(0.5, 5, 0.2), 4),
        "DE_L3_shallow_dense": round(seis.liquefaction_reduction_DE(0.8, 5, 0.4), 4),
        "DE_noliq_FL_ge_1": round(seis.liquefaction_reduction_DE(1.2, 5, 0.2), 4),
        "_note": "表8-1:第一級F_L≤1/3淺層鬆砂D_E=0(參數設零)/密砂1/6/深層1/3;第二級1/3~2/3;第三級2/3~1;F_L≥1不折減。深度10m為界,R_s>0.3密砂折減較輕。對齊公式卡_液狀化與基礎耐震"},
    "retrofit_shared_cracked": {
        "config": "R2/R4 同梁 b400/h01750/As1964(4D25)/C30 開裂換算斷面(α_Es=6.667,M_d1=200)",
        "x1_mm": round(retro_x1, 1), "Icr_e9mm4": round(retro_Icr / 1e9, 3),
        "eps_c1": round(retro_ec1, 6),
        "_note": "彈性開裂換算 x1=√(A1²+B1)−A1=191.3(非極限應力塊117);I_cr;ε_c1=M_d1·x1/(E_c·I_cr)。R2鋼板與R4增大截面同原梁→共用同一x1(納入引擎逼出R2算例x1修正)"},
    "retrofit_R1_CFRP": {
        "config": "R1 碳纖維CFRP抗彎(JTG/T J22) 碳布2層×寬300 t_f0.167/E_f2.4e5/ε_fu0.0155",
        "eps_f_allow": round(retro_epsf, 4), "km1": round(retro.cfrp_km1(2, 2.4e5, 0.167), 4),
        "xi_fb": round(retro_R1.xi_fb, 3), "x_mm": round(retro_R1.x, 1),
        "case2": retro_R1.case2, "Mu_kNm": round(retro_R1.Mu_kNm, 1),
        "_note": "[ε_f]=min(κ_m·ε_fu,⅔ε_fu,0.007)=0.007絕對上限控制;κ_m=min(κ_m1=0.813,κ_m2=0.85)≤0.9;ξ_fb=0.249;x=147.9≤ξ_fb·h→案②;M_u(式6-42)=539.4(+20%)。對齊算例_碳纖維CFRP抗彎補強設計"},
    "retrofit_R2_plate": {
        "config": "R2 外貼鋼板抗彎(JTG/T J22) 鋼板200×4 Q345 f_sp305/E_sp2.06e5",
        "x_mm": round(retro_R2.x, 1), "plate_yields": retro_R2.plate_yields,
        "sigma_sp_MPa": round(retro_R2.sigma_sp, 1), "Mu_kNm": round(retro_R2.Mu_kNm, 1),
        "l_p_mm": round(retro.plate_dev_length(305, 800, 2.5, 200), 1),
        "_note": "ε_sp需求1846≫f_sp→鋼板降伏取305;軸力平衡x=161.6;M_u(式6-26)=609(+36%,近40%上限);粘貼延伸l_p(式6-37)=788。對齊算例_外貼鋼板抗彎補強設計"},
    "cantilever_X_C2": _cantilever_C2(),
    "allowables_TW_audit": {
        "config": "稽核線批4/6（2026-09-24 NLM 7d947294）台灣明文容許值；f'c=40、f'ci=32、L=40000mm",
        "comp_service_TW_permanent": round(allowables.comp_service_TW_permanent(40.0), 4),
        "comp_service_TW_live_half": round(allowables.comp_service_TW_live_half(40.0), 4),
        "tension_precomp_bonded_general": round(allowables.tension_precompressed_TW(40.0), 4),
        "tension_precomp_bonded_segmental": round(allowables.tension_precompressed_TW(40.0, segmental=True), 4),
        "tension_precomp_corrosive_general": round(allowables.tension_precompressed_TW(40.0, corrosive=True), 4),
        "tension_precomp_corrosive_segmental": allowables.tension_precompressed_TW(40.0, segmental=True, corrosive=True),
        "tension_precomp_unbonded": allowables.tension_precompressed_TW(40.0, bonded=False),
        "transfer_tension_TW_unbonded": round(allowables.transfer_tension_TW(32.0), 4),
        "transfer_tension_TW_bonded": round(allowables.transfer_tension_TW(32.0, bonded=True), 4),
        "transfer_tension_TW_segmental_unbonded": allowables.transfer_tension_TW(32.0, segmental=True),
        "transfer_tension_current_AASHTO_nocap": round(allowables.transfer_tension(32.0), 4),
        "defl_limit": {k: round(allowables.deflection_limit_TW(40000.0, k), 4)
                       for k in allowables.TW_DEFLECTION_DENOM},
        "long_term_Ig": allowables.long_term_factor_TW(),
        "long_term_Ie_r05": round(allowables.long_term_factor_TW(False, 0.5), 4),
        "long_term_Ie_r15_floor": round(allowables.long_term_factor_TW(False, 1.5), 4),
        "_note": "🔴 comp_service_TW_permanent 原 −0.45f'c（AASHTO 5.9.4.2.1 值）→ 台灣 §8.15.2 2.(2) 明文 −0.40f'c，"
                 "更正時無呼叫者故既有 golden 不動。transfer_tension（現行 0.25√f'ci 無上限）與 "
                 "transfer_tension_TW（min(0.25√f'ci, 14 kgf/cm²=1.3729)）並列，f'ci=32 時差 3.0%，"
                 "是否將上限併入現行式待裁示。長期因素下限 1.6（§7.1.22 7.(4)b）A's/As>1.17 起控制。"},
    "retrofit_R4_enlargement": {
        "config": "R4 增大截面抗彎(JTG/T J22) 底加100mm(h800→900)+新筋2D25 h02860",
        "h0_mm": round((1964*750+982*860)/2946, 1), "x_mm": round(retro_R4.x, 1),
        "added_bar_yields": retro_R4.added_bar_yields, "sigma_s2_MPa": round(retro_R4.sigma_s2, 1),
        "Mu_kNm": round(retro_R4.Mu_kNm, 1),
        "_note": "新增筋ε_s2需求≫f_sd2→降伏取330;軸力平衡x=176.1;h0=A_s1+A_s2合力點786.7;M_u(式6-2)=679(+52%,超鋼板40%上限故選增大截面)。對齊算例_增大斷面補強設計"},
}

if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "golden_answers.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(golden, f, ensure_ascii=False, indent=2)
    print(f"已產生 {out}")
