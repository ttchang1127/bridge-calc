"""重新產生 golden_answers.json（Python 引擎與 JS 前端的共用驗證源）。

CI 在每次 push 時執行此檔，確保 golden 永遠來自最新的 bridgecalc。
執行：python3 make_golden.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bridgecalc import (Section, Tendon, compute_losses, combinations,
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
                        continuous_prestress)
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


def _continuous_pier():
    bot, top = _cont_case_groups()
    g = [bot, top]
    r = continuous_prestress([40, 40], g)
    M2p, M2m, M215 = r.X[1], r.M2_at(20), r.M2_at(15)
    e15 = bot.segs[0].e(15)
    e_eff = -(23700 * 80 + 12557 * 646) / 36257
    st_p, sb_p = pier_service_stress(36257e3, sec, e_eff, -41080 + M2p)
    _, sb_p0 = pier_service_stress(36257e3, sec, e_eff, -41080)
    st_15, sb_15 = pier_service_stress(23700e3, sec, e15, 14670 + 12000 + M215)
    _, sb_150 = pier_service_stress(23700e3, sec, e15, 14670 + 12000)
    Mu = 75337 - M2p
    ft = flexural_strength_T(11292, 1860, 40, 1400, 200, 700, 1975, Mu)
    tp = cont_tendon_segs([40, 40], -80, 950, -600, 0.12, 0.40)
    ra = continuous_prestress([40, 40], [TendonGroup(23700, tp.segs), top])
    _, sb_a15 = pier_service_stress(23700e3, sec, tp.e_at(15), 14670 + 12000 + ra.M2_at(15))
    return {
        "M1_x15_kNm": round(primary_moment_at(g, 15)), "M1_mid_kNm": round(primary_moment_at(g, 20)),
        "M1_pier_kNm": round(primary_moment_at(g, 40)),
        "M2_x15_kNm": round(M215), "M2_mid_kNm": round(M2m), "M2_pier_kNm": round(M2p),
        "M2_pier_bot_only_kNm": round(continuous_prestress([40, 40], [bot]).X[1]),
        "M2_pier_top_only_kNm": round(continuous_prestress([40, 40], [top]).X[1]),
        "int_M1_m_kNm2": round(r.b[0]), "flex_F_m": round(r.F[0][0], 3),
        "e_pier_eff_mm": round(e_eff, 1),
        "pier_service_sigma_top_MPa": round(st_p, 2), "pier_service_sigma_bot_MPa": round(sb_p, 2),
        "pier_service_sigma_bot_noM2_MPa": round(sb_p0, 2),
        "x15_service_sigma_top_MPa": round(st_15, 2), "x15_service_sigma_bot_MPa": round(sb_15, 2),
        "x15_service_sigma_bot_noM2_MPa": round(sb_150, 2),
        "pier_Mu_kNm": round(Mu), "pier_c_mm": round(ft.c), "pier_flanged": ft.flanged,
        "pier_fps_MPa": round(ft.fps), "pier_Mn_kNm": round(ft.Mn), "pier_phiMn_kNm": round(ft.phiMn),
        "pier_CR": round(ft.CR, 2), "pier_inadequate": not ft.ok,
        "adj_M2_pier_kNm": round(ra.X[1]), "adj_x15_sigma_bot_MPa": round(sb_a15, 2),
        "_note": "2026-09-17 校正：頂板腱 e 900→646(y_t−75−50)、底板腱垂度 1030、M2 改力法→正彎矩(跨中不利/墩頂有利)。"
                 "B墩底緣含M2 −12.70✓(不計M2 −19.72 為假性超限)；x=15m 底緣 +0.25 拉(台灣零拉✗)；"
                 "中墩 Mu=75,337−M2、dp=1975 → CR 0.55✗；調整：底板腱墩頂 −600 分段拋物線 → −2.02✓",
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
Vu_HS20 = 1419 + 229 + 681 * taiwan_per_lane_shear(40) / 588   # 活載剪力按 HS20/HL-93 縮放 ≈ 2,069
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
