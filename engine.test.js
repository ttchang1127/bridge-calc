/* engine.test.js — 單一引擎回歸：對 golden_answers.json 一次驗全部四域（= Python bridgecalc）。
 * 由原四份 *-engine.test.js 收斂；各域檢核塊包 IIFE 隔離區域變數、共用 chk/pass/fail。
 * 執行：node engine.test.js  （CI 亦跑；任一項不符即非零退出）。
 */
const { BC, SE, CE, RF } = require('./engine.js');
const g = require('./golden_answers.json');
var pass = 0, fail = 0;
function fmt(v) { return typeof v === 'boolean' ? v : (+v.toPrecision(6)); }
function chk(name, got, exp, tol) {
  var ok = Math.abs(got - exp) <= tol; ok ? pass++ : fail++;
  console.log((ok ? '\u2713' : '\u2717') + ' ' + name + ': JS ' + fmt(got) + ' vs golden ' + exp + (ok ? '' : '  \u274c \u0394>' + tol));
}
function chkEq(name, got, exp) {
  var ok = got === exp; ok ? pass++ : fail++;
  console.log((ok ? '\u2713' : '\u2717') + ' ' + name + ': JS ' + got + ' vs golden ' + exp + (ok ? '' : '  \u274c'));
}

// ═══════════════ 箱梁核心＋次檢核＋連續梁（原 box-girder-engine.test.js）═══════════════
(function () {
  // 40m 參考輸入（台灣 HS20-44 / 8組×19股）
  var sec = BC.section(5.065e6, 3.287e12, 1329, 2100);
  var t = BC.tendon(8, 19, 1109);
  var M_LL = BC.laneLiveLoad(BC.taiwanPerLaneMoment(40), 2, 1.0);
  var L = BC.computeLosses(t, sec, 24800, 4000);
  var c = BC.combinations(24800, 4000, M_LL);
  var s = BC.stresses(L.Pe, sec, t.e, c.Service_I);
  var Vu = 1419 + 229 + 681 * BC.taiwanPerLaneShear(40) / 588;
  var sh = BC.shearWeb(L.Pe, sec, t.e, 40, 250, 1692, Vu * 1e3, 1692, 40000);
  var fx = BC.flexuralStrength(t, sec, 40, 8000, 250, 1880, c.Strength_I, L.Pe, t.e);
  var df = BC.deflection(40000, 29700, sec, 144, L.Pe, t.e, 56.7 * BC.taiwanPerLaneMoment(40) / ((1 + 0.33) * 2867 + 1860));

  chk('影響線 HS20卡車M', BC.taiwanTruckMoment(40), g.influence_simple_40m.truck_absmax_kNm, 5);
  chk('影響線 每車道M', BC.taiwanPerLaneMoment(40), g.influence_simple_40m.per_lane_M_LL_IM_kNm, 5);
  chk('影響線 每車道V', BC.taiwanPerLaneShear(40), g.live_load_TW_HS20_40m.per_lane_V_kN, 3);
  chk('L0 M_LL 2車道', M_LL, g.loads.M_LL_IM_2lane_kNm, 5);
  chk('L0 Strength I', c.Strength_I, g.loads.StrengthI_kNm, 30);
  chk('L0 Service I', c.Service_I, g.loads.ServiceI_kNm, 5);
  chk('B 損失%', L.loss_pct * 100, g.prestress.loss_pct, 0.3);
  chk('B Pe(kN)', L.Pe / 1e3, g.prestress.Pe_kN, 50);
  chk('B Pe_min', BC.PeMinZeroTension(sec, t.e, c.Service_I) / 1e3, g.prestress.Pe_min_kN, 50);
  chk('C1 底緣σ', s.sb, g.service.sigma_bot_MPa, 0.06);
  chk('C1 頂緣σ', s.st, g.service.sigma_top_MPa, 0.06);
  chk('D1 fpc', sh.fpc, g.shear_D1.fpc_MPa, 0.05);
  chk('D1 Vcw', sh.Vcw / 1e3, g.shear_D1.Vcw_kN, 15);
  chk('D1 σ1', sh.sigma1, g.shear_D1.sigma1_MPa, 0.1);
  chk('M1 Mn', fx.Mn, g.flexure_M1.Mn_kNm, 100);
  chk('M1 CR', fx.CR, g.flexure_M1.CR, 0.02);
  chk('C2C3 δ_LL', df.d_LL, g.deflection.delta_LL_mm, 0.5);
  chk('C2C3 預拱', df.camber, g.deflection.camber_mm, 1);

  // ── 次檢核（v2）──
  var tr = BC.torsionCheck(L.Pe, sec, 40, 23.1e6, 26200, 1900);
  chk('D2 Tcr', tr.Tcr, g.torsion_D2.Tcr_kNm, 100);
  chk('D3 墩面 φMn', BC.slabFlexure(150.3, 2534, 200, 40, 420).phiMn, g.transverse_D3.support_phiMn_kNm, 0.5);
  var R_LL = 290 * BC.taiwanPerLaneShear(40) / 588;
  var br = BC.bearingCheck(1440 + R_LL, 1440, R_LL, 40, 100, 550, 450, 10, 8);
  chk('E1 σ_TL', br.sigma_TL, g.bearing_E1.sigma_TL_MPa, 0.05);
  chk('E1 形狀係數', br.shape_S, g.bearing_E1.shape_S, 0.1);
  var ej = BC.expansionJoint(8.8, 12.6, 8.0, 20);
  chk('E2 最大開度', ej.g_max, g.expansion_E2.g_max_mm, 0.1);
  var an = BC.anchorageCheck(t.Pi / 1e3, 8, 260, 2100, 4);
  chk('F1 Pu', an.Pu, g.anchorage_F1.Pu_kN, 5);
  chk('F1 剝落筋', an.As_spall, g.anchorage_F1.As_spall_mm2, 5);
  chk('F1 螺旋 Pult', BC.spiralLocalBearing(an.Pu, 2919, 8.47, 104044, 50, 380).Pult, g.anchorage_F1.spiral_Pult_kN, 5);
  var fa = BC.fatigueCheck(sec, L.Pe, t.e, 28800, 3222, 40);
  chk('P1 Δσ_ps', fa.dsig_ps, g.fatigue_P1.dsig_ps_MPa, 0.2);
  chk('P1 σ_c', fa.sig_c_max, g.fatigue_P1.sig_c_max_MPa, 0.1);
  chk('P1 箍筋@250', BC.stirrupFatigue(565, 250, 402, 1692).dfsv, g.fatigue_P1.stirrup_250_MPa, 3);
  var tbands = BC.thermalBandsFromDims(11000, 250, 5800, 200, 350, 2, 2100);
  var tg = BC.selfEquilibratingStress(tbands, 3.287e12, 771, 2100, [['底', 2100, 0]], 29700, 1.08e-5);
  chk('T1 TL', tg.TL, g.temperature_integrated_T1.TL_C, 0.2);
  chk('T1 負梯度底σ_SE', tg.sigma_neg['底'], g.temperature_integrated_T1.sigSE_bot_neg_MPa, 0.03);
  chk('T1 Service(含預力)', BC.thermalServiceCheck(tg.sigma_neg['底'], s.sb, 0.5).total, g.temperature_integrated_T1.service_total_MPa, 0.05);

  // ── 連續梁中墩（次彎矩 M2 + T 斷面）──
  var cp = g.continuous_pier;
  chk('連續 M2 跨中', BC.secondaryMoment(8320, BC.primaryMoment([[23700, 0.950], [12557, -0.300]])), cp.M2_mid_kNm, 5);
  chk('連續 M2 B墩', BC.secondaryMoment(-10594, BC.primaryMoment([[23700, -0.080], [12557, 0.900]])), cp.M2_pier_kNm, 5);
  var ft = BC.flexuralStrengthT(11292, 1860, 40, 1400, 200, 700, 1950, 75337);
  chk('中墩 T斷面 c', ft.c, cp.pier_c_mm, 2);
  chk('中墩 fps', ft.fps, cp.pier_fps_MPa, 5);
  chk('中墩 Mn', ft.Mn, cp.pier_Mn_kNm, 100);
  chk('中墩 CR(不足)', ft.CR, cp.pier_CR, 0.02);
})();

// ═══════════════ 耐震 S1/S2/S3/S5（原 seismic-engine.test.js）═══════════════
(function () {
  // ── S1 落橋防止（§8.5）：40m跨 / 基面起H=10m / 第二類地盤 Le=50m ──
  var s1 = g.seismic_S1;
  chk('S1 min N_L', SE.minFalloffLength(40, 10, 0), s1.min_NL_cm, 0.05);
  chk('S1 min N_L(斜角30°)', SE.minFalloffLength(40, 10, 30), s1.min_NL_skew30_cm, 0.01);
  chk('S1 地盤相對變位u_G', SE.groundRelativeDisplacement('第二類', 5000, 1.2), s1.u_G_cm, 0.01);
  chk('S1 防落長需求L_N', SE.requiredFalloffLength(70, 15, 22.5), s1.L_N_required_cm, 0.05);
  chk('S1 裝置降伏F_y', SE.restrainerYieldStrength(1800), s1.restrainer_Fy_kN, 0.5);

  // ── S2 隔震 LRB 等效線性化迭代（第7章）：W=8000/Q_d=400/K_d=6kN/mm/剛性墩/S_II,1=0.60 ──
  var s2 = g.seismic_S2, iso = SE.isolationDesign(8000, 400, 6000, 0.60);
  chk('S2 設計位移D_d(mm)', iso.Dd * 1000, s2.D_d_mm, 0.1);
  chk('S2 有效週期T_e', iso.Te, s2.T_e_s, 0.005);
  chk('S2 有效勁度K_eff', iso.Keff, s2.K_eff_kN_per_m, 0.5);
  chk('S2 系統阻尼ξ_e(%)', iso.xiE * 100, s2.xi_e_pct, 0.02);
  chk('S2 阻尼修正B_1', iso.B1, s2.B1, 0.0005);
  chk('S2 設計剪力V_b', iso.Vb_secant, s2.V_b_kN, 0.5);
  chkEq('S2 迭代收斂次數', iso.iterations, s2.iterations);
  chk('S2 B_1(15%)表3-1內插', SE.dampingCorrectionB1(0.15), s2.B1_15pct_table3_1, 0.0005);

  // ── S3 橋墩韌性容量設計（§4.2/5.3）：圓D150/fc280/fyh2800/Pe800000kgf/Mn3000tfm ──
  var s3 = g.seismic_S3;
  chk('S3 超強彎矩M_p', SE.overstrengthMoment(3000), s3.M_p_tfm, 0.5);
  chk('S3 容量剪力V_u', SE.capacityShear(3900, 8), s3.V_u_tf, 0.5);
  chk('S3 螺箍體積比ρ_s(%)', SE.rhoSCircular(280, 2800, 17671, 15394, 800000) * 100, s3.rho_s_pct, 0.002);
  chk('S3 圍束間距上限', SE.confinementSpacingLimit(150, 3.6), s3.spacing_limit_cm, 0.05);
  chk('S3 塑鉸長ℓ0', SE.plasticHingeLength(150, 800), s3.plastic_hinge_L_cm, 0.5);

  // ── S5 液狀化土壤參數折減 D_E（§8.1 表8-1）──
  var s5 = g.seismic_S5;
  chk('S5 D_E 一級/淺/鬆(參數設零)', SE.liquefactionReductionDE(0.3, 5, 0.2), s5.DE_L1_shallow_loose, 1e-4);
  chk('S5 D_E 一級/淺/密', SE.liquefactionReductionDE(0.3, 5, 0.4), s5.DE_L1_shallow_dense, 1e-4);
  chk('S5 D_E 一級/深', SE.liquefactionReductionDE(0.3, 15, 0.2), s5.DE_L1_deep, 1e-4);
  chk('S5 D_E 二級/淺/鬆', SE.liquefactionReductionDE(0.5, 5, 0.2), s5.DE_L2_shallow_loose, 1e-4);
  chk('S5 D_E 三級/淺/密', SE.liquefactionReductionDE(0.8, 5, 0.4), s5.DE_L3_shallow_dense, 1e-4);
  chk('S5 D_E 不液化(F_L≥1)', SE.liquefactionReductionDE(1.2, 5, 0.2), s5.DE_noliq_FL_ge_1, 1e-4);
})();

// ═══════════════ 施工 H1-H6（原 construction-engine.test.js）═══════════════
(function () {
  // ── H1/H2 支架施拉應力歷程：40m參考橋 8組×19股、f'ci=32、基準斷面 ──
  var h = g.construction_stage_H1H2;
  var sec = BC.section(5.065e6, 3.287e12, 1329, 2100), e = 1109, Pi = 29700e3, fci = 32;
  chk('H1/H2 施拉容許拉0.25√f\'ci', CE.transferTensionLimit(fci), h.transfer_tension_limit_MPa, 0.01);

  var h8 = CE.batchedTransfer(Pi, 8, 8, sec, e, fci);          // 全 PT 一次張拉（過平衡）
  chk('H2 全PT 頂緣σ_t', h8.st, h.S2_fullPT_top_MPa, 0.01);
  chkEq('H2 全PT 頂緣判定(超限)', h8.top_ok, h.S2_fullPT_top_ok);
  chk('H2 全PT 底緣σ_b', h8.sb, h.S2_fullPT_bot_MPa, 0.01);

  var h4 = CE.batchedTransfer(Pi, 4, 8, sec, e, fci);          // 分批 4/8 組
  chk('H2 分批4組 頂緣σ_t', h4.st, h.S2_batch4_top_MPa, 0.01);
  chkEq('H2 分批4組 頂緣判定(通過)', h4.top_ok, h.S2_batch4_top_ok);
  chk('H2 分批4組 底緣σ_b', h4.sb, h.S2_batch4_bot_MPa, 0.01);

  var h3s = CE.stageStress(Pi, sec, e, 24800, fci);            // 脫架（自重活化）
  chk('H2 脫架 頂緣σ_t', h3s.st, h.S3_strike_top_MPa, 0.01);
  chk('H2 脫架 底緣σ_b', h3s.sb, h.S3_strike_bot_MPa, 0.01);

  // ── H3 平衡懸臂：80+80m 變深 h_pier4.5/h_mid2.2、8節塊/側 ──
  var c = g.cantilever_H3;
  var w_h3 = [643, 599, 550, 497, 439, 385, 353, 341];
  var a_h3 = [6.75, 11.25, 15.75, 20.25, 24.75, 29.25, 33.75, 38.25].map(function (x) { return x - 4; });
  chk('H3 h_mid(x=0)', CE.variableDepth(0, 4.5, 2.2, 40), c.h_mid_m, 0.01);
  chk('H3 h_pier(x=半跨)', CE.variableDepth(40, 4.5, 2.2, 40), c.h_pier_m, 0.01);
  chk('H3 h(x=20)', CE.variableDepth(20, 4.5, 2.2, 40), c.h_at_x20_m, 0.001);
  chk('H3 自重彎矩ΣG·arm', CE.cantileverMoment(w_h3, a_h3), c.M_selfweight_about_x4_kNm, 1);
  chk('H3 懸臂最大M(含掛籃)', CE.cantileverMoment(w_h3, a_h3, 800, 40.5), c.M_cant_max_kNm, 1);
  chk('H3 長期下撓δ(1+φ)', CE.longTermDeflection(147, 2.0), c.delta_long_term_mm, 1);

  // ── H4 推進 ILM：40m等跨、等深 h=2.2m、A=4.870m²、Zb=3.093e9mm³ ──
  var l = g.launching_H4;
  var Mpos = CE.launchingSpanMoment(120.0, 40.0);
  var Pc = CE.centricPrestressRequired(Mpos, 3.093e9, 4.870e6, 1.5);
  chk('H4 懸臂根部M⁻', CE.launchingCantileverMoment(120.0, 14.0), l.M_cantilever_neg_kNm, 1);
  chk('H4 跨中M⁺', Mpos, l.M_span_pos_kNm, 1);
  chk('H4 臨時置中預力Pc', Pc, l.Pc_required_kN, 1);
  chkEq('H4 置中鋼腱束數', CE.nTendons(Pc, 2510), l.n_tendons);
  chk('H4 底緣σ_b(餘壓)', CE.launchingBottomStress(45100, 4.870e6, Mpos, 3.093e9), l.sigma_bot_min_MPa, 0.01);
  chk('H4 頂推力F=μW', CE.jackingForce(0.10, 120 * 200), l.jacking_force_kN, 1);
  chk('H4 滑動支承支壓σ_ba', CE.bearingStress(7200, 960000), l.bearing_stress_widened_MPa, 0.01);

  // ── H5/H6 預鑄節塊與接縫：H5 Ac=4.20m² / H6 接縫 Ac=2.85m² ──
  var s = g.segmental_H5H6;
  var Vkey = CE.shearKeyDesignCapacity(350, 0.439);
  chk('H5 節塊自重W=AcLγ', CE.segmentWeight(4.20, 2.5, 25), s.H5_segment_weight_kN, 0.1);
  chk('H5 接縫最小預力0.21Ac', CE.jointMinPrestress(4.20e6), s.H5_joint_min_prestress_kN, 0.1);
  chk('H5 接縫壓(4腱)', CE.jointCompression(4 * 480, 4.20e6), s.H5_joint_compression_4tendon_MPa, 0.001);
  chk('H6 接縫最小壓', CE.jointMinPrestress(2.85e6), s.H6_joint_min_compression_kN, 0.1);
  chk('H6 剪力鍵設計承載V_fuk·ξ', Vkey, s.H6_shear_key_design_capacity_kN, 0.1);
  chk('H6 LS3驗核比', CE.shearKeyUtilization(2850, 20, Vkey), s.H6_LS3_utilization, 0.01);
  chk('H6 黏結PT比例(≥30%)', CE.bondedPtRatio(14000, 42500), s.H6_bonded_pt_ratio, 0.001);
})();

// ═══════════════ 補強 R1/R2/R4（原 retrofit-engine.test.js）═══════════════
(function () {
  // ── 共用：彈性開裂換算斷面（二次受力）α_Es = 2.0e5/3.0e4 = 6.667、M_d1 = 200 kN·m ──
  var sc = g.retrofit_shared_cracked, aEs = 6.667;
  var x1 = RF.crackedNaDepth(400, 750, 1964, aEs);
  var Icr = RF.crackedInertia(400, x1, 1964, 750, aEs);
  var ec1 = RF.initialConcreteStrain(200, x1, 3.0e4, Icr);
  chk('共用 開裂中性軸x1', x1, sc.x1_mm, 0.05);
  chk('共用 開裂慣性矩I_cr(e9)', Icr / 1e9, sc.Icr_e9mm4, 0.005);
  chk('共用 初始混凝土應變ε_c1', ec1, sc.eps_c1, 5e-7);

  // ── R1 碳纖維 CFRP 抗彎（式6-42/6-44）：碳布2層×寬300、t_f0.167、E_f2.4e5、ε_fu0.0155 ──
  var r1g = g.retrofit_R1_CFRP;
  var epsf = RF.cfrpAllowableStrain(2, 2.4e5, 0.167, 0.0155);
  var r1 = RF.cfrpMomentCapacity(400, 800, 750, 13.8, 0.0033, 1964, 330, 100.2, 2.4e5, epsf, 0.0003);
  chk('R1 允許拉應變[ε_f](0.007控制)', epsf, r1g.eps_f_allow, 1e-5);
  chk('R1 κ_m1', RF.cfrpKm1(2, 2.4e5, 0.167), r1g.km1, 5e-5);
  chk('R1 界限相對受壓區ξ_fb', r1.xi_fb, r1g.xi_fb, 0.0005);
  chk('R1 受壓區高x', r1.x, r1g.x_mm, 0.05);
  chkEq('R1 案②(x≤ξ_fb·h)', r1.case2, r1g.case2);
  chk('R1 補強後M_u', r1.Mu_kNm, r1g.Mu_kNm, 0.05);

  // ── R2 外貼鋼板抗彎（式6-26/6-35/6-37）：鋼板200×4 Q345、f_sp305、E_sp2.06e5 ──
  var r2g = g.retrofit_R2_plate;
  var r2 = RF.plateMomentCapacity(400, 800, 750, 13.8, 0.0033, 1964, 330, 800, 305, 2.06e5, x1, ec1);
  chk('R2 受壓區高x', r2.x, r2g.x_mm, 0.05);
  chkEq('R2 鋼板降伏', r2.plate_yields, r2g.plate_yields);
  chk('R2 鋼板應力σ_sp', r2.sigma_sp, r2g.sigma_sp_MPa, 0.05);
  chk('R2 補強後M_u', r2.Mu_kNm, r2g.Mu_kNm, 0.05);
  chk('R2 粘貼延伸長度l_p', RF.plateDevLength(305, 800, 2.5, 200), r2g.l_p_mm, 0.05);

  // ── R4 增大截面抗彎（式6-2/6-10）：底加100mm(h800→900) + 新筋2D25、h02=860 ──
  var r4g = g.retrofit_R4_enlargement;
  var h0c = (1964 * 750 + 982 * 860) / 2946;
  var r4 = RF.enlargementMomentCapacity(400, 860, h0c, 13.8, 0.0033, 1964, 330, 982, 330, 2.0e5, x1, ec1);
  chk('R4 加固後有效高h0', h0c, r4g.h0_mm, 0.05);
  chk('R4 受壓區高x', r4.x, r4g.x_mm, 0.05);
  chkEq('R4 新增筋降伏', r4.added_bar_yields, r4g.added_bar_yields);
  chk('R4 新增筋應力σ_s2', r4.sigma_s2, r4g.sigma_s2_MPa, 0.05);
  chk('R4 補強後M_u', r4.Mu_kNm, r4g.Mu_kNm, 0.05);
})();

console.log('\n' + pass + '/' + (pass + fail) + ' \u5c0d\u9f4a golden' +
  (fail ? ' \u274c ' + fail + ' \u9805\u4e0d\u7b26' : ' \u2713 JS \u5f15\u64ce = Python bridgecalc\uff08\u56db\u57df\u5168\u9a57\uff09'));
process.exit(fail ? 1 : 0);
