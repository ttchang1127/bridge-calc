/* box-girder-engine.js 回歸測試：JS 引擎對 golden_answers.json 驗證（= Python bridgecalc）。
 * 執行：node box-girder-engine.test.js   （CI 亦跑；不符即非零退出）
 */
var BC = require('./box-girder-engine.js');
var g = require('./golden_answers.json');
var pass = 0, fail = 0;
function chk(name, got, exp, tol) {
  var ok = Math.abs(got - exp) <= tol; ok ? pass++ : fail++;
  console.log((ok ? '✓' : '✗') + ' ' + name + ': JS ' + (+got.toFixed(2)) + ' vs golden ' + exp + (ok ? '' : '  ❌ Δ>' + tol));
}

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

console.log('\n' + pass + '/' + (pass + fail) + ' 對齊 golden' + (fail ? ' ❌ ' + fail + ' 項不符' : ' ✓ JS 引擎 = Python 引擎（核心 + 次檢核）'));
process.exit(fail ? 1 : 0);
