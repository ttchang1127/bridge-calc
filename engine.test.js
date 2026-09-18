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
  var Vu = BC.taiwanContShearAt([40], 1.692, 'R', 5.065 * 24.5, 20, 2).Vu_pos / 2;   // d_v 斷面實算÷2 腹板 ≈ 2,298
  var sh = BC.shearWeb(L.Pe, sec, t.e, 40, 250, 1692, Vu * 1e3, 1692, 40000);
  var fx = BC.flexuralStrength(t, sec, 40, 8000, 250, 1880, c.Strength_I, L.Pe, t.e);
  var df = BC.deflection(40000, 29700, sec, 144, L.Pe, t.e, 56.7 * BC.taiwanPerLaneMoment(40) / ((1 + 0.33) * 2867 + 1860));

  chk('影響線 HS20卡車M', BC.taiwanTruckMoment(40), g.influence_simple_40m.truck_absmax_kNm, 5);
  chk('影響線 每車道M', BC.taiwanPerLaneMoment(40), g.influence_simple_40m.per_lane_M_LL_IM_kNm, 5);
  chk('影響線 每車道V', BC.taiwanPerLaneShear(40), g.live_load_TW_HS20_40m.per_lane_V_kN, 3);
  chk('影響線 卡車支承V（雙向）', BC.taiwanTruckShear(40), g.live_load_TW_HS20_40m.truck_V_support_kN, 0.1);
  chk('影響線 30m 每車道V（卡車控制）', BC.taiwanPerLaneShear(30), g.live_load_TW_HS20_40m.per_lane_V_30m_kN, 0.2);
  chk('L0 M_LL 2車道', M_LL, g.loads.M_LL_IM_2lane_kNm, 5);
  chk('L0 Strength I', c.Strength_I, g.loads.StrengthI_kNm, 30);
  chk('L0 Service I', c.Service_I, g.loads.ServiceI_kNm, 5);
  chk('B 損失%', L.loss_pct * 100, g.prestress.loss_pct, 0.3);
  chk('B Pe(kN)', L.Pe / 1e3, g.prestress.Pe_kN, 50);
  chk('B Pe_min', BC.PeMinZeroTension(sec, t.e, c.Service_I) / 1e3, g.prestress.Pe_min_kN, 50);
  // 反解設計（design.py 移植）：設計＝驗算之逆，對 design_inverse golden
  var Pe_min = BC.PeMinZeroTension(sec, t.e, c.Service_I);
  var Pe21 = BC.computeLosses(BC.tendon(8, 21, 1109), sec, 24800, 4000).Pe;
  var w_DL = 8 * (24800 + 4000) * 1e6 / (40000 * 40000);
  chk('設計 Pe_min(反解種子)', Pe_min / 1e3, g.design_inverse.Pe_min_zero_tension_kN, 1);
  chk('設計 最小鋼腱組數', BC.minTendonGroups(Pe_min, 19, L.fpe), g.design_inverse.min_tendon_groups, 0);
  chk('設計 所需垂度(LBR0.985)', BC.requiredDrape(0.985, w_DL, 40000, Pe21), g.design_inverse.required_drape_LBR0985_mm, 2);
  chk('設計 最小 Sb', BC.minSectionModulusSb(L.Pe, sec.A, t.e, c.Service_I, 0) / 1e9, g.design_inverse.min_Sb_zero_tension_e9mm3, 0.01);
  // G1 管道實配（duct_layout 移植）：對三組 golden ＋ 形心守恆自檢
  (function () {
    chk('G1 曲率半徑 R(m)', BC.radiusOfCurvature(1109, 40000) / 1000, g.tendon_profile_G1.R_m, 0.5);
    // 摩擦損失沿長度分佈：對 friction_profile_G1 golden ＋ 與 tendon_profile_G1 交叉驗證
    (function () {
      var fg = g.friction_profile_G1, sg = BC.parabolicCurvSegs(40, 1109);
      chk('G1 摩擦 κ', sg[0][2], fg.kappa_rad_per_m, 1e-6);
      chk('G1 摩擦 α跨中', BC.frictionAngle(sg, 20, 40, false), fg.alpha_mid_rad, 1e-4);
      chk('G1 摩擦 α全長', BC.frictionAngle(sg, 40, 40, false), fg.alpha_full_rad, 1e-4);
      ['both', 'start', 'end', 'alt'].forEach(function (jk) {
        var p = BC.frictionProfile(40, sg, 0.25, 0.003, jk), e = fg[jk];
        chk('G1 摩擦 ' + jk + ' 起點', p.at_start, e.start, 5e-4);
        chk('G1 摩擦 ' + jk + ' 跨中', p.at_mid, e.mid, 5e-4);
        chk('G1 摩擦 ' + jk + ' 終點', p.at_end, e.end, 5e-4);
        chk('G1 摩擦 ' + jk + ' 平均', p.avg, e.avg, 5e-4);
        chk('G1 摩擦 ' + jk + ' 最大', p.max_ratio, e.max, 5e-4);
        chk('G1 摩擦 ' + jk + ' 最大位置', p.x_max, e.x_max_m, 0.1);
      });
      // 與 tendon_profile_G1 交叉：同一組 μ/K/a/L 由兩條路徑得同值
      chk('G1 摩擦 跨中≡dual_mid', BC.frictionAt(20, 40, sg, 0.25, 0.003, 'both'),
          g.tendon_profile_G1.friction_dual_mid, 1e-3);
      chk('G1 摩擦 遠端≡single_end', BC.frictionAt(40, 40, sg, 0.25, 0.003, 'start'),
          g.tendon_profile_G1.friction_single_end, 1e-3);
      // 接進 computeLosses：不傳＝原式；傳了就改變 Pe
      var L0 = BC.computeLosses(t, sec, 24800, 4000);
      var L1 = BC.computeLosses(t, sec, 24800, 4000, { fricRatio: BC.frictionAt(20, 40, sg, 0.25, 0.003, 'both') });
      var L2 = BC.computeLosses(t, sec, 24800, 4000, { fricRatio: BC.frictionAt(40, 40, sg, 0.25, 0.003, 'start') });
      chk('G1 摩擦 接入(跨中雙端≡預設)', L1.Pe / 1e3, L0.Pe / 1e3, 3);
      chk('G1 摩擦 接入(單端遠端 ΔPe kN)', (L0.Pe - L2.Pe) / 1e3,
          t.fpj * (fg.start.end - fg.both.mid) * t.Aps / 1e3, 5);
    })();
    // 逐腱合成：對 tendon_forces_G1 golden
    (function () {
      var tg = g.tendon_forces_G1, sg = BC.parabolicCurvSegs(40, 1109);
      function tf(rows, mode) {
        var jk = BC.assignJack(rows.length, mode), ts = rows.map(function (r, i) {
          return { no: 'T' + (i + 1), y: r[0], x: r[1], jack: jk[i] }; });
        return BC.tendonForces(ts, 10, 40, sg, sec.yb, 1395, 140 * 19, 180.0, 0.25, 0.003);
      }
      var A = [[-80, 0], [120, 0], [320, 0], [520, 0]], B = [[150, -85], [150, 85], [290, -85], [290, 85]];
      var a = tf(A, 'alt'), b = tf(B, 'alt'), st = tf(A, 'start');
      chk('G1 逐腱 單列4層 Pe', a.Pe_total / 1e3, tg.single_col_4row.Pe_kN, 0.2);
      chk('G1 逐腱 單列4層 e_eff', a.e_eff, tg.single_col_4row.e_eff_mm, 0.05);
      chk('G1 逐腱 單列4層 Δe', a.de, tg.single_col_4row.de_mm, 0.05);
      chk('G1 逐腱 單列4層 x_eff', a.x_eff, tg.single_col_4row.x_eff_mm, 0.05);
      chk('G1 逐腱 2列2層 Δe', b.de, tg.two_col_2row.de_mm, 0.05);
      chk('G1 逐腱 2列2層 x_eff', b.x_eff, tg.two_col_2row.x_eff_mm, 0.05);
      chk('G1 逐腱 全自起點 Pe', st.Pe_total / 1e3, tg.single_col_start.Pe_kN, 0.2);
      chk('G1 逐腱 全自起點 Δe', st.de, tg.single_col_start.de_mm, 1e-6);
      chk('G1 逐腱 總力≡平均損失率', a.Pe_total / 1e3, a.Pe_avg_ratio / 1e3, 0.002);
      chkEq('G1 逐腱 assignJack', BC.assignJack(4, 'alt').join(','), 'start,end,start,end');
    })();
    chk('G1 淨距需求 tw', BC.ductSpacingRequired(100, 25), 40, 0);
    chk('G1 淨距需求 od', BC.ductSpacingRequired(100, 25, 'od'), 100, 0);
    var dg = g.duct_layout_G1;
    var d = BC.ductLayout(8, 2, sec.yb - 1109, { od: 100, webT: 350, cover: 40, sv: 40, dAgg: 25, yb: sec.yb, h: 2100 });
    chkEq('G1 實配 列數', d.nCol, dg.n_col);
    chkEq('G1 實配 層數', d.nRow, dg.n_row);
    chk('G1 實配 淨間距需求', d.sReq, dg.s_req_mm, 0.1);
    chk('G1 實配 腹板可用寬', d.webAvail, dg.web_avail_mm, 0.1);
    chk('G1 實配 層距', d.pitchV, dg.pitch_v_mm, 0.1);
    chk('G1 實配 水平淨距', d.sH, dg.s_h_mm, 0.1);
    chk('G1 實配 最底層 y', d.yBot, dg.y_bot_mm, 0.1);
    chk('G1 實配 底層管外緣', d.coverBot, dg.cover_bot_mm, 0.1);
    chkEq('G1 實配 保護層 OK', d.coverOk, dg.cover_ok);
    chkEq('G1 實配 可排下', d.fits, dg.fits);
    chk('G1 實配 e_max', d.eMax, dg.e_max_mm, 0.5);
    d.rows.forEach(function (r, i) { chk('G1 實配 第' + (i + 1) + '層 y', r.y, dg.rows_y_mm[i], 0.1); });

    var sg = g.duct_layout_G1_strict;                     // 保守側：同配置改排不下
    var ds = BC.ductLayout(8, 2, sec.yb - 1109, { od: 100, webT: 350, cover: 75, sv: 100, dAgg: 25, yb: sec.yb, h: 2100, rule: 'od' });
    chkEq('G1 保守 列數', ds.nCol, sg.n_col);
    chkEq('G1 保守 層數', ds.nRow, sg.n_row);
    chk('G1 保守 底層管外緣', ds.coverBot, sg.cover_bot_mm, 0.1);
    chkEq('G1 保守 可排下', ds.fits, sg.fits);
    chk('G1 保守 e_max', ds.eMax, sg.e_max_mm, 0.5);
    chk('G1 單點簡化 e_max', ds.eMaxPoint, g.duct_layout_G1_strict.e_max_mm + 300, 0.5);

    var pg = g.duct_layout_G1_pier;                       // 墩頂：頂緣控制
    var dp = BC.ductLayout(8, 2, sec.yb + 900, { od: 100, webT: 350, cover: 40, sv: 40, dAgg: 25, yb: sec.yb, h: 2100 });
    chkEq('G1 墩頂 穿出頂緣', dp.topOk, pg.top_ok);
    chk('G1 墩頂 頂層 y', dp.yTop, pg.y_top_mm, 0.1);
    chk('G1 墩頂 e_pier 下限', sec.yb - 2100 + 40 + 50 + (dp.yTop - dp.yCgs), pg.e_pier_limit_mm, 0.5);

    // 形心守恆：3D 畫多腱不得改變引擎在用的偏心 e（含餘數層的不對稱填充）
    [[8, 2, 350, 100], [5, 1, 500, 100], [12, 3, 600, 60], [7, 2, 450, 50]].forEach(function (cs) {
      var q = BC.ductLayout(cs[0], cs[1], 260, { od: 90, webT: cs[2], cover: 50, sv: cs[3], dAgg: 20, rule: 'od' }), tot = 0, cnt = 0;
      q.rows.forEach(function (r) { tot += r.y * r.count; cnt += r.count; });
      chk('G1 實配 形心守恆 n=' + cs[0] + '/web' + cs[1], tot / cnt, 260, 1e-6);
    });
  })();
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

  // ── 連續梁（次彎矩 M2 力法 + 服務性 + 中墩 T 斷面）──
  var cp = g.continuous_pier;
  var cb = -(950 + 80) / 400, ctp = (300 + 646) / 225;
  var gBot = { P: 23700, segs: [BC.parabolaSeg(0, 40, 20, 950, cb), BC.parabolaSeg(40, 80, 60, 950, cb)] };
  var gTop = { P: 12557, segs: [BC.parabolaSeg(25, 40, 40, -646, ctp), BC.parabolaSeg(40, 55, 40, -646, ctp)] };
  var fmc = BC.continuousPrestress([40, 40], [gBot, gTop]);
  chk('連續 M1 跨中', BC.primaryMomentAt([gBot, gTop], 20), cp.M1_mid_kNm, 1);
  chk('連續 M1 B墩', BC.primaryMomentAt([gBot, gTop], 40), cp.M1_pier_kNm, 1);
  chk('連續 M2 x=15（力法）', fmc.M2At(15), cp.M2_x15_kNm, 1);
  chk('連續 M2 跨中（力法）', fmc.M2At(20), cp.M2_mid_kNm, 1);
  chk('連續 M2 B墩（力法）', fmc.X[1], cp.M2_pier_kNm, 1);
  chk('連續 ∫M1·m', fmc.b[0], cp.int_M1_m_kNm2, 1);
  chk('連續 M2 底板腱單獨', BC.continuousPrestress([40, 40], [gBot]).X[1], cp.M2_pier_bot_only_kNm, 1);
  chk('連續 M2 頂板腱單獨', BC.continuousPrestress([40, 40], [gTop]).X[1], cp.M2_pier_top_only_kNm, 1);
  var envC = BC.taiwanContEnvelope([40, 40], 5.065 * 24.5, 20, 2, 40), pierC = envC[40];
  chk('中墩 M_DC（解析）', pierC.M_dc, cp.pier_M_dc_kNm, 1);
  chk('中墩 M_LL−（HS20-44 2車道）', pierC.M_ll_neg, cp.pier_M_ll_neg_kNm, 1);
  chk('中墩 衝擊（相鄰兩跨平均）', pierC.I_neg, cp.pier_impact, 1e-4);
  var sbP = BC.stresses(36257e3, BC.section(5.065e6, 3.287e12, 1329, 2100), cp.e_pier_eff_mm, pierC.Ms_neg + fmc.X[1]).sb;
  chk('中墩 σb（含M2）', sbP, cp.pier_service_sigma_bot_MPa, 0.01);
  var ft = BC.flexuralStrengthT(4 * 2660, 1860, 40, 1400, 200, 700, 1975, -(pierC.Mu_neg + fmc.X[1]));
  chk('中墩 Mu(含M2)', -(pierC.Mu_neg + fmc.X[1]), cp.pier_Mu_kNm, 1);
  chk('中墩 T斷面 c', ft.c, cp.pier_c_mm, 1);
  chk('中墩 fps', ft.fps, cp.pier_fps_MPa, 1);
  chk('中墩 Mn', ft.Mn, cp.pier_Mn_kNm, 1);
  chk('中墩 CR(含M2)', ft.CR, cp.pier_CR, 0.01);
  // 連續梁解析影響線＋台灣 HS20-44 包絡
  var ce = g.cont_envelope_taiwan;
  chk('連續IL 墩頂 p=20', BC.contMomentIL([40, 40], 40, 20), ce.il_pier_p20, 1e-6);
  chk('連續IL x=15 p=15', BC.contMomentIL([40, 40], 15, 15), ce.il_x15_p15, 1e-6);
  chk('連續DL 墩頂', BC.contDLMoment([40, 40], 130.4, 40), ce.dl_pier_w130_4, 1e-3);
  chk('連續DL x=15', BC.contDLMoment([40, 40], 130.4, 15), ce.dl_x15_w130_4, 1e-3);
  var lvP = BC.taiwanContLiveMoment([40, 40], 40), lv15 = BC.taiwanContLiveMoment([40, 40], 15);
  chk('墩頂 卡車負彎矩', lvP.truck_neg, ce.pier_truck_neg, 0.01);
  chk('墩頂 車道負彎矩（2集中）', lvP.lane_neg, ce.pier_lane_neg, 0.01);
  chk('墩頂 衝擊', lvP.I_neg, ce.pier_I_neg, 1e-6);
  chk('x15 卡車正', lv15.truck_pos, ce.x15_truck_pos, 0.01);
  chk('x15 車道正', lv15.lane_pos, ce.x15_lane_pos, 0.01);
  chk('x15 卡車負', lv15.truck_neg, ce.x15_truck_neg, 0.01);
  chk('x15 車道負', lv15.lane_neg, ce.x15_lane_neg, 0.01);
  chk('單跨退化 車道', BC.taiwanContLiveMoment([40], 20).lane_pos, ce.single_span_lane_mid, 0.01);
  chk('單跨退化 卡車', BC.taiwanContLiveMoment([40], 20).truck_pos, ce.single_span_truck_mid, 0.01);
  var lv3 = BC.taiwanContLiveMoment([30, 40, 30], 30);
  chk('三跨 墩1 車道負', lv3.lane_neg, ce.three_span_pier1_lane_neg, 0.01);
  chk('三跨 墩1 卡車負', lv3.truck_neg, ce.three_span_pier1_truck_neg, 0.01);
  chk('三跨 墩1 衝擊(L=35)', lv3.I_neg, ce.three_span_pier1_I_neg, 1e-6);
  var env20 = BC.taiwanContEnvelope([40, 40], 124.0925, 20, 2, 20);
  chk('包絡 min Mu−', Math.min.apply(null, env20.map(function (r) { return r.Mu_neg; })), ce.env_min_Mu_neg_kNm, 0.1);
  chk('包絡 max Mu+', Math.max.apply(null, env20.map(function (r) { return r.Mu_pos; })), ce.env_max_Mu_pos_kNm, 0.1);
  chk('多車道折減 3', BC.taiwanLaneReduction(3), ce.lane_reduction_3, 1e-9);
  chk('多車道折減 4', BC.taiwanLaneReduction(4), ce.lane_reduction_4, 1e-9);
  // 雙系統配束：墩頂局部腱＋頂板構造檢核
  var pcg = g.pier_cap_tendon, pc = BC.pierCapTendonSegs([40, 40], 15, 300, -646);
  chk('頂板腱 c', pc.segs[0].c, pcg.c_mm_per_m2, 1e-6);
  chk('頂板腱 R_min', pc.Rmin, pcg.R_min_mm, 1);
  chk('頂板腱 錨點1', pc.anchors[0], pcg.anchors_m[0], 1e-9);
  chk('頂板腱 M2 單獨', BC.continuousPrestress([40, 40], [{ P: 12557, segs: pc.segs }]).X[1], pcg.M2_top_only_kNm, 1);
  var c75 = BC.topSlabTendonCheck(2100, 1329, 250, -646, 4, 5100, 100, 75), c40 = BC.topSlabTendonCheck(2100, 1329, 250, -646, 4, 5100, 100, 40);
  var cbad = BC.topSlabTendonCheck(2100, 1329, 250, -700, 12, 5100, 100, 40, 25, 'od');
  chk('頂板 e_hi(75)', c75.e_hi, pcg.slab75_e_hi, 1e-9); chk('頂板 e_lo(75)', c75.e_lo, pcg.slab75_e_lo, 1e-9);
  chkEq('頂板 ok(75)', c75.ok, pcg.slab75_ok);
  chk('頂板 e_hi(40)', c40.e_hi, pcg.slab40_e_hi, 1e-9); chk('頂板 e_lo(40)', c40.e_lo, pcg.slab40_e_lo, 1e-9);
  chk('頂板 水平淨距', c40.s_clear, pcg.slab40_s_clear, 1e-9);
  chkEq('頂板 超出(in_slab)', cbad.in_slab, pcg.bad_in_slab);
  chk('頂板 超出 頂保護層', cbad.top_cover, pcg.bad_top_cover, 1e-9);
  chk('頂板 12束淨距', cbad.s_clear, pcg.bad_s_clear, 0.01);
  var tpD = BC.contTendonSegs([40, 40], 0, 1109, -600);
  chk('雙系統 合成 M2', BC.continuousPrestress([40, 40], [{ P: 23724, segs: tpD.segs }, { P: 12557, segs: pc.segs }]).X[1], pcg.dual_M2_pier_kNm, 0.1);
  chk('三跨 頂板腱長受限', BC.pierCapTendonSegs([30, 40, 30], 20, 300, -600).lengths[0], pcg.three_span_lengths[0], 1e-9);
  // 變斷面（中墩底板加厚）：數值柔度、階梯 EI 解析、M2 與包絡
  var vs = g.variable_section_haunch, dimsV = [11000, 250, 5800, 200, 350, 2, 2100];
  var prV = BC.haunchProfile.apply(null, [[40, 40]].concat(dimsV, [400, 8])), pr0V = BC.haunchProfile.apply(null, [[40, 40]].concat(dimsV, [200, 8]));
  var refV = BC.sectionFromDims.apply(null, dimsV), s40V = prV.sectionAt(40);
  chk('變斷面 基準 I', refV.I / 1e12, vs.ref_I_e12, 1e-5);
  chk('變斷面 墩頂 A', s40V.A, vs.pier_A, 1);
  chk('變斷面 墩頂 I', s40V.I / 1e12, vs.pier_I_e12, 1e-5);
  chk('變斷面 墩頂 ȳb', s40V.yb, vs.pier_yb, 1e-3);
  chk('變斷面 Δȳ', prV.dyb(40), vs.dyb_pier, 1e-3);
  chk('變斷面 I_rel', prV.Irel(40), vs.I_rel_pier, 1e-6);
  chk('變斷面 底板厚 x=36', prV.botTAt(36), vs.bot_t_x36, 1e-9);
  var tpV = BC.contTendonSegs([40, 40], 0, 1109, -600), gV = [{ P: 23724, segs: tpV.segs }];
  chk('數值柔度＝閉合（M2）', BC.continuousPrestress([40, 40], gV, 8, pr0V).X[1], vs.numeric0_M2_pier, 0.01);
  chk('數值柔度＝閉合（包絡）', BC.taiwanContEnvelope([40, 40], 124.0925, 20, 2, 20, 0.25, 400, pr0V)[20].Mu_neg, vs.numeric0_env_pier_Mu_neg, 0.1);
  var stV = BC.contFlex([40, 40], function (x) { return Math.abs(x - 40) < 8 ? 2 : 1; }, [32, 48]);
  chk('階梯 EI 均布（解析）', stV.supportMomentsDist(function () { return 10; })[1], vs.step_EI_dl_pier_w10, 1e-3);
  chk('階梯 EI 點載重', stV.supportMomentsPoint(20)[1], vs.step_EI_point_p20, 1e-6);
  chk('加厚 M2（全長腱）', BC.continuousPrestress([40, 40], gV, 8, prV).X[1], vs.haunch_M2_pier_default, 0.1);
  chk('加厚 M2（雙系統）', BC.continuousPrestress([40, 40], [gBot, gTop], 8, prV).X[1], vs.haunch_M2_pier_dual, 0.1);
  var ehV = BC.taiwanContEnvelope([40, 40], 124.0925, 20, 2, 20, 0.25, 400, prV);
  chk('加厚 墩頂 M_DC', ehV[20].M_dc, vs.haunch_env_pier_M_dc, 0.1);
  chk('加厚 墩頂 M_LL−', ehV[20].M_ll_neg, vs.haunch_env_pier_M_ll_neg, 0.1);
  chk('加厚 墩頂 Mu−', ehV[20].Mu_neg, vs.haunch_env_pier_Mu_neg, 0.1);
  chk('加厚 x16 Mu+', ehV[8].Mu_pos, vs.haunch_env_x16_Mu_pos, 0.1);
  // 連續梁剪力：影響線、DL、HS20-44、次剪力、算例墩左 d_v 斷面抗剪
  var csg = g.cont_shear_taiwan;
  chk('剪力IL x10R p30', BC.contShearIL([40, 40], 10, 30, 'R'), csg.il_x10R_p30, 1e-8);
  chk('剪力IL 墩左 p20', BC.contShearIL([40, 40], 40, 20, 'L'), csg.il_pierL_p20, 1e-8);
  chk('剪力IL 墩右 p60', BC.contShearIL([40, 40], 40, 60, 'R'), csg.il_pierR_p60, 1e-8);
  chk('DL 剪力 墩左', BC.contDLShear([40, 40], 10, 40, 'L'), csg.dl_pierL_w10, 1e-6);
  chk('DL 剪力 端', BC.contDLShear([40, 40], 10, 0, 'R'), csg.dl_end_w10, 1e-6);
  var lvLs = BC.taiwanContLiveShear([40, 40], 40, 'L'), lvSs = BC.taiwanContLiveShear([40], 0, 'R'), lv3s = BC.taiwanContLiveShear([30, 40, 30], 30, 'R');
  chk('墩左 卡車剪力', lvLs.truck_neg, csg.pierL_truck_neg, 0.01);
  chk('墩左 車道剪力', lvLs.lane_neg, csg.pierL_lane_neg, 0.01);
  chk('單跨退化 卡車剪力', lvSs.truck_pos, csg.single_truck, 0.01);
  chk('單跨退化 車道剪力', lvSs.lane_pos, csg.single_lane, 0.01);
  chk('三跨 墩1右 車道剪力', lv3s.lane_pos, csg.three_pier1R_lane_pos, 0.01);
  chk('三跨 墩1右 衝擊', lv3s.I, csg.three_pier1R_I, 1e-6);
  var rowS = BC.taiwanContShearAt([40, 40], csg.case_x_m, 'L', 5.065 * 24.5, 20, 2);
  chk('算例墩左 V_DC', rowS.V_dc, csg.case_V_dc, 0.01);
  chk('算例墩左 V_LL−', rowS.V_ll_neg, csg.case_V_ll_neg, 0.01);
  chk('算例墩左 Vu−', rowS.Vu_neg, csg.case_Vu_neg, 0.01);
  var V2s = BC.secondaryShear(fmc, [40, 40], csg.case_x_m, 'L'), VuS = BC.designShearWithV2(rowS.Vu_pos, rowS.Vu_neg, V2s);
  chk('算例墩左 V2', V2s, csg.case_V2, 0.01);
  chk('算例墩左 設計Vu', VuS, csg.case_Vu_design, 0.01);
  var shS = BC.shearWebAt(36257e3, BC.section(5.065e6, 3.287e12, 1329, 2100), csg.case_P_slope_kN / 36257, 40, 250, 1512, VuS / 2 * 1e3, 2);
  chk('算例墩左 Vp/腹板', shS.Vp / 1e3, csg.case_Vp_per_web_kN, 0.01);
  chk('算例墩左 σ1', shS.sigma1, csg.case_sigma1, 0.001);
  chk('算例墩左 Vcw', shS.Vcw / 1e3, csg.case_Vcw_kN, 0.01);
  chk('算例墩左 Av/s 需', shS.Av_s_req, csg.case_Av_s_req, 1e-4);
  chk('算例墩左 φVn @250', BC.phiVn(shS.Vcw, 397.4 / 250, 1512) / 1e3, csg.case_phiVn_D16x2_s250_kN, 0.1);
  chk('算例墩左 φVn @200', BC.phiVn(shS.Vcw, 397.4 / 200, 1512) / 1e3, csg.case_phiVn_D16x2_s200_kN, 0.1);
  // 頂板腱逐點損失（摩擦＋滑移＋ES＋長期）
  var ptl = g.pier_cap_tendon_losses, pcL = BC.pierCapTendonSegs([40, 40], 15, 300, -646);
  var f40L = BC.pierCapTendonForce(pcL, 40, 1395, 10640, ptl.other), f25L = BC.pierCapTendonForce(pcL, 25, 1395, 10640, ptl.other);
  chk('頂板腱 墩頂摩擦', f40L.friction, ptl.pier_friction, 0.01);
  chk('頂板腱 墩頂 f_pe', f40L.fpe, ptl.pier_fpe, 0.01);
  chk('頂板腱 墩頂 P', f40L.P, ptl.pier_P_kN, 0.1);
  chk('頂板腱 錨碇滑移', f25L.slip, ptl.anchor_slip, 0.01);
  chk('頂板腱 錨碇 f_pe', f25L.fpe, ptl.anchor_fpe, 0.01);
  chk('滑移（L_set 超過半長）', BC.anchorSlipLoss(0, 15, 5, 6).dsigma, ptl.slip_capped_dsigma, 1e-3);
  var topL = { P: function (x) { return BC.pierCapTendonForce(pcL, x, 1395, 10640, ptl.other).P; }, segs: pcL.segs };
  chk('頂板腱精算 M2 墩頂', BC.continuousPrestress([40, 40], [gBot, topL]).X[1], ptl.M2_pier_kNm, 0.1);
  // 全長腱錨具滑移
  var tsg = g.tendon_slip_G1, fpjS = 0.75 * 1860, seg40S = BC.parabolicCurvSegs(40, 1109);
  chk('滑移 簡支40 L_set', Math.sqrt(6 * 195000 / 1000 / (fpjS * BC.frictionAt(20, 40, seg40S, 0.25, 0.003, 'start') / 20)), tsg.simple40_L_set_m, 1e-3);
  chk('滑移 簡支40 錨端', BC.tendonSlipLoss(0, 40, seg40S, fpjS, 'both', 0.25, 0.003, 6), tsg.simple40_slip_anchor, 0.01);
  chk('滑移 簡支40 跨中（＝0）', BC.tendonSlipLoss(20, 40, seg40S, fpjS, 'both', 0.25, 0.003, 6), tsg.simple40_slip_mid, 1e-9);
  var tpS = BC.contTendonSegs([40, 40], 0, 1109, -600), csS = tpS.segs.map(function (gg) { return [gg.x1, gg.x2, 2 * Math.abs(gg.c) / 1000]; });
  chk('滑移 連續80 起端張拉 x=1.5', BC.tendonSlipLoss(1.5, 80, csS, fpjS, 'start', 0.25, 0.003, 6), tsg.cont80_slip_start_jack_x1_5, 0.01);
  chk('滑移 連續80 終端張拉 x=1.5', BC.tendonSlipLoss(1.5, 80, csS, fpjS, 'end', 0.25, 0.003, 6), tsg.cont80_slip_end_jack_x1_5, 1e-9);
  chk('滑移 預設 0（不計）', BC.tendonSlipLoss(1.5, 80, csS, fpjS, 'start', 0.25, 0.003, 0), 0, 1e-12);
  // 中間錨碇段
  var mag = g.mid_anchor_G1, tpM = BC.contTendonSegs([40, 40], 0, 1109, -600);
  var csM = tpM.segs.map(function (gg) { return [gg.x1, gg.x2, 2 * Math.abs(gg.c) / 1000]; }), fpjM = 0.75 * 1860;
  chk('一段 80m 最大損失%', BC.segmentedFrictionProfile([0, 80], csM, fpjM, 0.25, 0.003, 'both', 0).max_ratio * 100, mag.one_seg_max_pct, 0.01);
  chk('中間錨碇 2 段 最大%', BC.segmentedFrictionProfile([0, 40, 80], csM, fpjM, 0.25, 0.003, 'both', 6).max_ratio * 100, mag.two_seg_slip6_max_pct, 0.01);
  chk('中間錨碇 2 段（純摩擦）', BC.segmentedFrictionProfile([0, 40, 80], csM, fpjM, 0.25, 0.003, 'both', 0).max_ratio * 100, mag.two_seg_fric_only_max_pct, 0.01);
  chk('中間錨碇 3 段 最大%', BC.segmentedFrictionProfile([0, 80 / 3, 160 / 3, 80], csM, fpjM, 0.25, 0.003, 'both', 6).max_ratio * 100, mag.three_seg_slip6_max_pct, 0.01);
  var fmM = BC.segmentedTendonForce(40, [0, 40, 80], csM, fpjM, 8 * 19 * 140, 0, 0.25, 0.003, 'both', 6);
  chk('中間錨碇 墩頂 f_pe', fmM.fpe, mag.pier_fpe_two_seg, 0.1);
  chk('中間錨碇 墩頂滑移', fmM.slip, mag.pier_slip_two_seg, 0.1);
  // 沿長度的預力 P(x)（loss_profile）
  var lpg = g.loss_profile_G1, secLP = BC.section(5.065e6, 3.287e12, 1329, 2100),
      tenLP = BC.tendon(8, 19, 1109), Llp = 40, aLp = 870,
      segLP = BC.parabolicCurvSegs(Llp, aLp),
      MDlp = BC.udlMoment(5.065e6 / 1e6 * 24.5, Llp),
      eLp = BC.parabolicE(tenLP.e, tenLP.e - aLp, Llp);
  var lp0 = BC.lossProfile(tenLP, secLP, Llp, eLp, MDlp, segLP, { jack: 'both', slip: 0, n: 40 }),
      lp1 = BC.lossProfile(tenLP, secLP, Llp, eLp, MDlp, segLP, { jack: 'both', slip: 6, n: 40 }),
      p0lp = lp1.at(0), pmlp = lp1.at(20);
  chk('P(x) 端部摩擦%', p0lp.fric / tenLP.fpj * 100, lpg.end_fric_pct, 0.01);
  chk('P(x) 端部滑移%', p0lp.slip / tenLP.fpj * 100, lpg.end_slip_pct, 0.01);
  chk('P(x) 端部損失%', p0lp.loss_pct * 100, lpg.end_loss_pct, 0.01);
  chk('P(x) 端部 Pe', p0lp.Pe / 1e3, lpg.end_Pe_kN, 0.1);
  chk('P(x) 跨中摩擦%', pmlp.fric / tenLP.fpj * 100, lpg.mid_fric_pct, 0.01);
  chk('P(x) 跨中滑移%（＝0）', pmlp.slip / tenLP.fpj * 100, lpg.mid_slip_pct, 1e-9);
  chk('P(x) 跨中損失%', pmlp.loss_pct * 100, lpg.mid_loss_pct, 0.01);
  chk('P(x) 跨中 Pe', pmlp.Pe / 1e3, lpg.mid_Pe_kN, 0.1);
  chk('P(x) 無滑移 x_Pemin（跨中）', lp0.x_Pemin, lpg.x_Pemin_noslip_m, 1e-9);
  chk('P(x) 有滑移 x_Pemin（端部）', lp1.x_Pemin, lpg.x_Pemin_slip_m, 1e-9);
  chk('P(x) Pe_min 無滑移', lp0.Pe_min / 1e3, lpg.Pe_min_noslip_kN, 0.1);
  chk('P(x) Pe_min 有滑移', lp1.Pe_min / 1e3, lpg.Pe_min_slip_kN, 0.1);
  chk('P(x) 全長平均 Pe', lp1.Pe_avg / 1e3, lpg.Pe_avg_slip_kN, 0.1);
  var singLP = BC.computeLosses(tenLP, secLP, MDlp(Llp / 2), 0,
                                { fricRatio: BC.frictionAt(Llp / 2, Llp, segLP, 0.25, 0.003, 'both') });
  chk('P(x) 單點式跨中 Pe（對照）', singLP.Pe / 1e3, lpg.single_point_mid_Pe_kN, 0.1);
  chk('P(x) 單點式在 x=1m 高估%',
      (singLP.Pe - lp1.at(1).Pe) / lp1.at(1).Pe * 100, lpg.single_point_overestimate_at_1m_pct, 0.01);
  // 中間錨碇齒塊錨碇區
  var blg = g.blister_B1, bld = BC.blisterDesign(1764, 3000, 3000,
      { a_plate: 200, b_plate: 200, L_b: 400, W_b: 300, D_b: 250, A_bearing: 60000,
        fci: 35, f_cb: 8, A_cb: 40000, alpha_deg: 5, mu: 1, fy: 420, fsd: 360, fc: 40,
        straight_have: 400 });
  chk('齒塊 錨板承壓 f_b', bld.bearing.f_b, blg.f_b_MPa, 0.05);
  chk('齒塊 承壓容許值', bld.bearing.f_b_allow, blg.f_b_allow_MPa, 0.01);
  chkEq('齒塊 承壓通過（需螺旋筋）', bld.bearing.ok, blg.bearing_ok);
  chk('齒塊 螺旋筋需提升倍數', bld.bearing.spiral_factor_req, blg.spiral_factor_req, 0.01);
  chk('齒塊 爆裂力 F_burst', bld.burst.F_burst, blg.F_burst_kN, 0.05);
  chk('齒塊 爆裂合力位置 d_burst', bld.burst.d_burst, blg.d_burst_mm, 0.05);
  chk('齒塊 爆裂配筋 As', bld.burst.As_burst, blg.As_burst_mm2, 0.1);
  chk('齒塊 Tie-back 0.25Ps', bld.tie.T_req, blg.T_tieback_kN, 0.05);
  chk('齒塊 既有預壓抵扣', bld.tie.C_precomp, blg.C_precomp_kN, 0.05);
  chk('齒塊 Tie-back 容許應力', bld.tie.fs_allow, blg.fs_allow_MPa, 0.05);
  chk('齒塊 Tie-back As（抵扣）', bld.tie.As, blg.As_tie_mm2, 0.1);
  chk('齒塊 Tie-back As（保守）', bld.tie.As_conservative, blg.As_tie_conservative_mm2, 0.1);
  chk('齒塊 剝裂力', bld.spall.F_spall, blg.F_spall_kN, 0.01);
  chk('齒塊 剝裂配筋', bld.spall.As_spall, blg.As_spall_mm2, 0.1);
  chk('齒塊 介面剪力 V_int', bld.face.V_int, blg.V_interface_kN, 0.05);
  chk('齒塊 介面剪力摩擦配筋', bld.face.As_vf, blg.As_vf_mm2, 0.1);
  chk('齒塊 介面剪應力 τ', bld.face.tau, blg.tau_interface_MPa, 0.001);
  chk('齒塊 介面剪應力上限', bld.face.tau_cap, blg.tau_cap_MPa, 0.01);
  chk('齒塊 介面剪力上限 V_cap', bld.face.V_cap, blg.V_cap_kN, 0.1);
  chkEq('齒塊 介面面積足夠（🔴算例漏檢）', bld.face.area_ok, blg.interface_area_ok);
  chk('齒塊 所需最小介面面積', bld.face.A_req, blg.A_interface_req_mm2, 1);
  chkEq('齒塊 控制配筋項', bld.governing, blg.governing);
  chk('齒塊 四類配筋合計', bld.As_total_conservative, blg.As_total_conservative_mm2, 0.1);
  chkEq('齒塊 幾何檢核', bld.geom.ok, blg.geom_ok);
  // 施工階段體系轉換的潛變重分配
  var stg = g.staging_redist_S1, spSt = [40, 40], wSt = 5.065 * 24.5, dphiSt = 1.7;
  var dSt = BC.spanBySpanDeadLoad(spSt, wSt, 40);
  var rSt = BC.creepRedistribution(dSt.xs, dSt.M_I, dSt.M_II, dphiSt, 0.8, 'trost');
  var rStD = BC.creepRedistribution(dSt.xs, dSt.M_I, dSt.M_II, dphiSt, 0.8, 'dischinger');
  var midSt = rSt.at(20), pierSt = rSt.at(40);
  chk('重分配 λ Trost', rSt.factor.lam, stg.lam_trost, 1e-4);
  chk('重分配 λ Dischinger', rStD.factor.lam, stg.lam_dischinger, 1e-4);
  chkEq('重分配 Trost Δφ=10 逸出範圍', !BC.redistributionFactor(10, 0.8, 'trost').valid,
        stg.lam_trost_invalid_at_dphi10);
  chk('逐跨 跨中 M_I（簡支）', midSt.M_I, stg.mid_M_I_kNm, 0.1);
  chk('逐跨 跨中 M_II（連續）', midSt.M_II, stg.mid_M_II_kNm, 0.1);
  chk('逐跨 跨中 M(∞)', midSt.M_inf, stg.mid_M_inf_kNm, 0.1);
  chk('逐跨 跨中／連續比', midSt.M_inf / midSt.M_II, stg.mid_ratio_to_cont, 1e-3);
  chk('逐跨 墩頂 M_I（＝0）', pierSt.M_I, stg.pier_M_I_kNm, 1e-6);
  chk('逐跨 墩頂 M_II（連續）', pierSt.M_II, stg.pier_M_II_kNm, 0.1);
  chk('逐跨 墩頂 M(∞)', pierSt.M_inf, stg.pier_M_inf_kNm, 0.1);
  chk('逐跨 墩頂／連續比', pierSt.M_inf / pierSt.M_II, stg.pier_ratio_to_cont, 1e-3);
  chkEq('🔑 束制彎矩支承間線性', BC.redistributionIsLinear(rSt, spSt), stg.restraint_linear);
  var cpSt = BC.continuousPrestress(spSt, [{ P: 23725, segs: BC.simpleSpanTendonSegs(spSt, 1109) }]),
      xpSt = [], iSt;
  for (iSt = 0; iSt <= 80; iSt++) xpSt.push(iSt);
  var m2cSt = xpSt.map(function (x) { return cpSt.M2At(x); });
  var rmSt = BC.prestressM2Redistribution(spSt, xpSt, m2cSt, dphiSt);
  chk('M₂ 墩頂（連續體系張拉）', rmSt.M2_pier_cont, stg.M2_pier_cont_kNm, 0.1);
  chk('M₂ 墩頂（先簡支後連續）', rmSt.M2_pier_inf, stg.M2_pier_inf_kNm, 0.1);
  chk('墩頂合計 一次成形', pierSt.M_II + rmSt.M2_pier_cont, stg.pier_total_cont_kNm, 0.1);
  chk('墩頂合計 逐跨施工', pierSt.M_inf + rmSt.M2_pier_inf, stg.pier_total_sbs_kNm, 0.1);
  chk('M₂ 簡支線形＝閉合解 P·a', rmSt.M2_pier_cont, stg.M2_pier_closed_form_Pa_kNm, 0.1);
  // 體系轉換套進連續梁包絡
  var s2 = g.staged_envelope_S2, lamS2 = 1.7 / (1 + 0.8 * 1.7),
      envS2 = BC.taiwanContEnvelope(spSt, wSt, 20, 2, 20, 0.25, 400, null),
      xsS2 = envS2.map(function (r) { return r.x; }), ipS2 = 0, k2;
  for (k2 = 1; k2 < xsS2.length; k2++) if (Math.abs(xsS2[k2] - 40) < Math.abs(xsS2[ipS2] - 40)) ipS2 = k2;
  var M2S2 = xsS2.map(function (x) { return cpSt.M2At(x); }),
      aS2 = BC.stagedEnvelope(spSt, envS2, wSt, lamS2),
      bS2 = BC.stagedEnvelope(spSt, envS2, wSt, lamS2, { M2: M2S2, psAtSimple: true });
  var imo = envS2.reduce(function (a, b) { return b.Mu_pos > a.Mu_pos ? b : a; }),
      ima = aS2.reduce(function (a, b) { return b.Mu_pos > a.Mu_pos ? b : a; });
  chk('包絡 一次成形 跨中 Mu⁺', imo.Mu_pos, s2.mono_mid_Mu_pos_kNm, 0.1);
  chk('包絡 一次成形 中墩 Mu⁻', envS2[ipS2].Mu_neg, s2.mono_pier_Mu_neg_kNm, 0.1);
  chk('包絡 逐跨 跨中 Mu⁺', ima.Mu_pos, s2.sbs_mid_Mu_pos_kNm, 0.1);
  chk('包絡 逐跨 跨中位置', ima.x, s2.sbs_mid_x_m, 1e-3);
  chkEq('包絡 逐跨 跨中控制狀態', ima.gov_pos, s2.sbs_mid_gov);
  chk('包絡 逐跨 中墩 Mu⁻', aS2[ipS2].Mu_neg, s2.sbs_pier_Mu_neg_kNm, 0.1);
  chkEq('包絡 逐跨 中墩控制狀態', aS2[ipS2].gov_neg, s2.sbs_pier_gov_neg);
  chk('包絡 簡支張拉 中墩 Mu⁻', bS2[ipS2].Mu_neg, s2.pss_pier_Mu_neg_kNm, 0.1);
  chk('包絡 簡支張拉 中墩 Mu⁺（正彎矩接頭）', bS2[ipS2].Mu_pos, s2.pss_pier_Mu_pos_kNm, 0.1);
  chkEq('包絡 簡支張拉 中墩 Mu⁺ 控制狀態', bS2[ipS2].gov_pos, s2.pss_pier_gov_pos);
  var timSt = BC.timingSensitivity(2.0, [['7d', 0.25], ['28d', 0.55], ['180d', 1.35]],
                                   0, -wSt * 1600 / 8);
  chk('時序 λ(7天合龍)', timSt[0].lam, stg.lam_7d, 1e-4);
  chk('時序 λ(180天合龍)', timSt[2].lam, stg.lam_180d, 1e-4);
  chk('時序 墩頂 M(7天)', timSt[0].M_pier, stg.M_pier_7d_kNm, 0.1);
  chk('時序 墩頂 M(180天)', timSt[2].M_pier, stg.M_pier_180d_kNm, 0.1);
  // 套管尺寸（§8.25.4／表 8.3／PTI 4.3）
  var dsg = g.duct_size_T83, dsA = BC.ductSizeCheck(19, 100), dsB = BC.ductSizeCheck(19, 90, null, 15.2, 140, 'pti_pull');
  chk('套管 面積比 19×15.2 φ100', dsA.ratio, dsg.ratio_19x152_id100, 1e-4);
  chkEq('套管 §8.25.4 面積≥2倍', dsA.area_ok, dsg.area_ok_tw);
  chk('套管 表8.3 最大內徑', dsA.id_max_tw, dsg.id_max_tw_19x152_mm, 1e-9);
  chkEq('套管 未給外徑＝把內徑當外徑', dsA.od_gt_id, dsg.od_gt_id_when_od_unknown);
  chk('套管 面積比 φ90', dsB.ratio, dsg.ratio_id90, 1e-4);
  chkEq('套管 φ90 PTI 拉入法 2.5 倍', dsB.area_ok, dsg.area_ok_id90_pti_pull);
  chkEq('套管 φ105 超表8.3', BC.ductSizeCheck(19, 105).id_ok, dsg.id_ok_id105);
  // 包絡恆載因數依有利／不利取 max／min
  var gmg = g.cont_envelope_gamma_min, envG = BC.taiwanContEnvelope([40, 40], 5.065 * 24.5, 20, 2, 20, 0.25, 400, null),
      atG = function (x) { return envG.reduce(function (a, b) { return Math.abs(b.x - x) < Math.abs(a.x - x) ? b : a; }); };
  chk('γ_min x=28 Mu⁻（負彎矩區延伸）', atG(28).Mu_neg, gmg.x28_Mu_neg_kNm, 0.1);
  chk('γ_min x=32 Mu⁺（墩旁正彎矩）', atG(32).Mu_pos, gmg.x32_Mu_pos_kNm, 0.1);
  chk('γ_min 墩頂 Mu⁺', atG(40).Mu_pos, gmg.pier_Mu_pos_kNm, 0.1);
  chk('γ_min 跨中 Mu⁻', atG(20).Mu_neg, gmg.x20_Mu_neg_kNm, 0.1);
  // 連續橫隔梁正彎矩接頭
  var pmg = g.pos_moment_conn_S4, pmA = BC.positiveMomentConnection(262.9, 3.287e12, 1329, 40, null, null, 2100 - 80),
      pmB = BC.positiveMomentConnection(262.9, 3.287e12, 1329, 40, 90, null, 2100 - 80);
  chk('正彎矩接頭 f_r', pmA.fr, pmg.fr_MPa, 1e-4);
  chk('正彎矩接頭 M_cr', pmA.Mcr, pmg.Mcr_kNm, 0.1);
  chkEq('正彎矩接頭 控制', pmA.governs, pmg.gov);
  chk('正彎矩接頭 需求', pmA.M_req, pmg.M_req_kNm, 0.1);
  chk('正彎矩接頭 As 估算', pmA.As_est, pmg.As_est_mm2, 1);
  chkEq('正彎矩接頭 90 天簡化 控制', pmB.governs, pmg.gov_90d);
  chk('正彎矩接頭 90 天簡化 需求', pmB.M_req, pmg.M_req_90d_kNm, 0.1);
  // 台灣第十二章耐久性保護層
  var dct = g.durability_cover_T12;
  chk('保護層 一般50年 II w/c0.45 表列', BC.durabilityCover('general', 50, 'II', 0.45).table, dct.gen50_II_wc045_table, 1e-9);
  chk('保護層 一般50年 II 需求（含§8.25.1底線）', BC.durabilityCover('general', 50, 'II', 0.45).required, dct.gen50_II_wc045_req, 1e-9);
  chk('保護層 一般100年 II', BC.durabilityCover('general', 100, 'II', 0.45).required, dct.gen100_II_wc045_req, 1e-9);
  chk('保護層 一般50年 I（箱內）', BC.durabilityCover('general', 50, 'I', 0.45).table, dct.gen50_I_wc045_table, 1e-9);
  chk('保護層 鹽害50年 嚴重 腹版', BC.durabilityCover('salt', 50, null, 0.40, '梁腹版外露面', '嚴重').required, dct.salt50_web_severe_req, 1e-9);
  chk('保護層 鹽害100年 極嚴重 腹版', BC.durabilityCover('salt', 100, null, 0.40, '梁腹版外露面', '極嚴重').required, dct.salt100_web_extreme_req, 1e-9);
  chk('保護層 橋面版頂層 鹽害中度50年', BC.durabilityCover('salt', 50, null, 0.45, '橋面版頂層筋', '中度').required, dct.deck_top_salt50_moderate, 1e-9);
  // duct_layout 側向額外餘裕
  var dsm = g.duct_side_margin, dlo = function (m) { return BC.ductLayout(8, 2, 1329 - 1109, { yb: 1329, h: 2100, sideMargin: m }); };
  chk('側向餘裕0 cover_side', dlo(0).coverSide, dsm.m0_cover_side, 1e-6);
  chk('側向餘裕0 s_h', dlo(0).sH, dsm.m0_s_h, 1e-6);
  chk('側向餘裕15 cover_side', dlo(15).coverSide, dsm.m15_cover_side, 1e-6);
  chk('側向餘裕15 s_h', dlo(15).sH, dsm.m15_s_h, 1e-6);
  chkEq('側向餘裕15 列數', dlo(15).nCol, dsm.m15_ncol);
  chkEq('側向餘裕30 列數（退單列）', dlo(30).nCol, dsm.m30_ncol);
  chk('側向餘裕30 e_max', dlo(30).eMax, dsm.m30_e_max, 0.1);
  chkEq('側向餘裕30 排得下', dlo(30).fits, dsm.m30_fits);
  // AASHTO 潛變係數＋嚴格 AAEM
  var cag = g.creep_aashto_A1, vsA = BC.boxVolumeSurface(11000, 250, 5800, 200, 350, 2, 2100),
      tsA = BC.timingSensitivityAashto(7, [28, 90, 180], 0, -5.065 * 24.5 * 1600 / 8, 75, vsA, 32);
  chk('AASHTO V/S（箱內周長 50%）', vsA, cag.VS_mm, 0.01);
  chk('AASHTO ψ(∞,7)', BC.aashtoCreep(Infinity, 7, 75, vsA, 32).psi, cag.psi_inf_t0_7, 1e-4);
  chk('28天合龍 Δφ', tsA[0].dphi, cag.t28_dphi, 1e-4);
  chk('28天合龍 φ(∞,t1)', tsA[0].phi_r, cag.t28_phi_r, 1e-4);
  chk('28天合龍 λ 嚴格', tsA[0].lam_exact, cag.t28_lam_exact, 1e-4);
  chk('28天合龍 λ 近似', tsA[0].lam_approx, cag.t28_lam_approx, 1e-4);
  chk('28天合龍 墩頂 M', tsA[0].M_pier, cag.t28_M_pier_kNm, 0.1);
  chk('90天合龍 λ 嚴格', tsA[1].lam_exact, cag.t90_lam_exact, 1e-4);
  chk('90天合龍 墩頂 M', tsA[1].M_pier, cag.t90_M_pier_kNm, 0.1);
  chk('180天合龍 λ 嚴格', tsA[2].lam_exact, cag.t180_lam_exact, 1e-4);
  // §8.25.3 套管捆紮
  var dbg = g.duct_bundle_825_3, db350 = BC.ductLayoutBundled(8, 2, 1329 - 1109, { webT: 350, yb: 1329, h: 2100 }),
      db300 = BC.ductLayoutBundled(8, 2, 1329 - 1109, { webT: 300, yb: 1329, h: 2100 });
  chk('捆紮 350 不捆 e_max', db350.cands.none.eMax, dbg.w350_none_e_max, 0.1);
  chkEq('捆紮 350 最佳型式', db350.best.bundle, dbg.w350_best);
  chk('捆紮 350 最佳 e_max', db350.best.eMax, dbg.w350_best_e_max, 0.1);
  chk('捆紮 300 不捆 e_max', db300.cands.none.eMax, dbg.w300_none_e_max, 0.1);
  chkEq('捆紮 300 不捆 可行', db300.cands.none.fits, dbg.w300_none_fits);
  chkEq('捆紮 300 最佳型式', db300.best.bundle, dbg.w300_best);
  chk('捆紮 300 最佳 e_max', db300.best.eMax, dbg.w300_best_e_max, 0.1);
  chkEq('捆紮 300 最佳 可行', db300.best.fits, dbg.w300_best_fits);
  // 簡支 d_v 斷面設計剪力（analyzer ⑤ 自動帶入 Vu）
  var ssd = g.simple_shear_dv, rSS = BC.taiwanContShearAt([40], ssd.x_m, 'R', 5.065 * 24.5, 20, 2);
  chk('簡支 d_v V_DC', rSS.V_dc, ssd.V_dc, 1e-3);
  chk('簡支 d_v V_LL+IM', rSS.V_ll_pos, ssd.V_ll, 1e-3);
  chk('簡支 d_v 衝擊（至較遠支點）', rSS.I, ssd.I, 1e-6);
  chk('簡支 d_v Vu/腹板', rSS.Vu_pos / 2, ssd.Vu_per_web, 0.01);
  // 連續梁全長剪力掃描＋箍筋分區
  var css = g.cont_shear_scan;
  var scan = BC.contShearDesignScan([40, 40], BC.groupsPrestressAt([gBot, gTop]), fmc, 5.065 * 24.5, 20, 2, 2100, 1329, 40, 250, 2, 397.4, 1.0,
                                    [24.95, 25.05, 54.95, 55.05], null, BC.section(5.065e6, 3.287e12, 1329, 2100));
  chkEq('剪力掃描 斷面數', scan.rows.length, css.n_rows);
  var worstS = scan.rows.reduce(function (a, b) { return b.Av_s_req > a.Av_s_req ? b : a; });
  chk('剪力掃描 最不利 x', worstS.x, css.worst_x, 1e-3);
  chk('剪力掃描 最不利 Av/s', worstS.Av_s_req, css.worst_Av_s, 1e-4);
  chkEq('剪力掃描 最不利間距', worstS.s_pick, css.worst_s);
  var aL = scan.rows.find(function (r) { return Math.abs(r.x - 24.95) < 1e-6; }), aR = scan.rows.find(function (r) { return Math.abs(r.x - 25.05) < 1e-6; });
  chk('錨碇左 Vp', aL.Vp_web, css.anchorL_Vp, 0.1); chk('錨碇右 Vp', aR.Vp_web, css.anchorR_Vp, 0.1);
  chkEq('錨碇左 間距', aL.s_pick, css.anchorL_s); chkEq('錨碇右 間距', aR.s_pick, css.anchorR_s);
  chkEq('分區數', scan.zones.length, css.zones.length);
  chkEq('分區內容', JSON.stringify(scan.zones.map(function (z) { return [Math.round(z[0] * 1000) / 1000, Math.round(z[1] * 1000) / 1000, z[2]]; })), JSON.stringify(css.zones));
  // 體系轉換的恆載剪力＋箍筋分區
  var sss = g.staged_shear_S3, lamSS = 1.7 / (1 + 0.8 * 1.7);
  var scanS = BC.contShearDesignScan([40, 40], BC.groupsPrestressAt([gBot, gTop]), fmc, 5.065 * 24.5, 20, 2, 2100, 1329, 40, 250, 2, 397.4, 1.0,
                                     [24.95, 25.05, 54.95, 55.05], null, BC.section(5.065e6, 3.287e12, 1329, 2100), null, null, lamSS);
  chk('逐跨 端支承 d_v Vu/腹板（一次成形）', scan.rows[0].Vu_web, sss.end_dv_Vu_web_mono, 0.1);
  chk('逐跨 端支承 d_v Vu/腹板（逐跨）', scanS.rows[0].Vu_web, sss.end_dv_Vu_web_sbs, 0.1);
  chkEq('逐跨 端部箍筋間距', scanS.zones[0][2], sss.end_zone_s_sbs);
  chk('逐跨 端部分區長', scanS.zones[0][1], sss.end_zone_to_sbs, 1e-3);
  chkEq('逐跨 分區內容', JSON.stringify(scanS.zones.map(function (z) { return [Math.round(z[0] * 1000) / 1000, Math.round(z[1] * 1000) / 1000, z[2]]; })), JSON.stringify(sss.zones_sbs));
  chk('§8.20.3 上限（減半）', BC.stirrupMaxSpacingTW(1e6, 40, 250, 1512, 2100)[0], css.max_spacing_halved, 1e-9);
  chk('§8.20.3 上限（一般）', BC.stirrupMaxSpacingTW(5e5, 40, 250, 1512, 2100)[0], css.max_spacing_normal, 1e-9);
  // HS20-44 中後軸距 4.25～9.15 掃描
  var ras = g.rear_axle_scan, rs = BC.taiwanRearSpacings(), l15 = BC.taiwanContLiveMoment([15, 15], 15), s40 = BC.taiwanContLiveMoment([40], 20);
  chkEq('軸距候選數', rs.length, ras.spacings_n);
  chk('軸距上限', rs[rs.length - 1], ras.spacings_last, 1e-9);
  chk('15+15 墩頂卡車（掃描）', l15.truck_neg, ras.p15_truck_neg_scan, 0.001);
  chk('15+15 控制軸距', l15.V_neg, ras.p15_V_neg, 1e-9);
  chk('簡支40 x20 卡車', s40.truck_pos, ras.simple40_x20_truck, 0.001);
  chk('簡支40 控制軸距', s40.V_pos, ras.simple40_x20_V, 1e-9);
  // 分析器預設全長連續腱（分段拋物線）＋三跨
  var ctd = g.cont_tendon_force_default, tpd = BC.contTendonSegs([40, 40], 0, 1109, -600);
  var fmd = BC.continuousPrestress([40, 40], [{ P: 23724, segs: tpd.segs }]);
  chkEq('線形段數', tpd.segs.length, ctd.n_segs);
  chk('反曲點1', tpd.infl[0], ctd.infl_m[0], 0.001);
  chk('反曲點2', tpd.infl[1], ctd.infl_m[1], 0.001);
  chk('線形 R_min', tpd.Rmin, ctd.R_min_mm, 1);
  chk('線形 e(10)', tpd.eAt(10), ctd.e_at_10, 0.01);
  chk('線形 e(36)', tpd.eAt(36), ctd.e_at_36, 0.01);
  chk('預設 M2 B墩（力法）', fmd.X[1], ctd.M2_pier_kNm, 0.1);
  chk('預設 M2 x=20', fmd.M2At(20), ctd.M2_x20_kNm, 0.1);
  var tp3 = BC.contTendonSegs([30, 40, 30], 0, 900, -500);
  var fm3 = BC.continuousPrestress([30, 40, 30], [{ P: 20000, segs: tp3.segs }]);
  chk('三跨 M2 墩1', fm3.X[1], ctd.three_span_X_kNm[1], 0.1);
  chk('三跨 M2 墩2', fm3.X[2], ctd.three_span_X_kNm[2], 0.1);
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
  chk('H3 懸臂最大M(含掛籃)', CE.cantileverMoment(w_h3, a_h3, 800, 36.5), c.M_cant_max_kNm, 1);
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
