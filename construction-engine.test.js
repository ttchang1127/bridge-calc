/* construction-engine.js 回歸測試：JS 引擎對 golden_answers.json 驗證
 * （= Python bridgecalc.construction + launching + segmental）。
 * 執行：node construction-engine.test.js  （CI 亦跑；不符即非零退出）
 *
 * ⚠️ H1/H2 須用「基準斷面」BC.section(5.065e6, 3.287e12, 1329, 2100)——與 golden/
 *    box-girder-engine.test.js 同源。若改用 sectionFromDims（I=3.2867e12）會差 ~0.01 MPa。
 */
var BC = require('./box-girder-engine.js');   // 供 CE.stageStress 取用 BC.stresses
global.BC = BC;
var CE = require('./construction-engine.js');
var g = require('./golden_answers.json');
var pass = 0, fail = 0;
function fmt(v) { return typeof v === 'boolean' ? v : (+v.toPrecision(6)); }
function chk(name, got, exp, tol) {
  var ok = Math.abs(got - exp) <= tol; ok ? pass++ : fail++;
  console.log((ok ? '✓' : '✗') + ' ' + name + ': JS ' + fmt(got) + ' vs golden ' + exp + (ok ? '' : '  ❌ Δ>' + tol));
}
function chkEq(name, got, exp) {
  var ok = got === exp; ok ? pass++ : fail++;
  console.log((ok ? '✓' : '✗') + ' ' + name + ': JS ' + got + ' vs golden ' + exp + (ok ? '' : '  ❌'));
}

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

console.log('\n' + pass + '/' + (pass + fail) + ' 對齊 golden' +
  (fail ? ' ❌ ' + fail + ' 項不符' : ' ✓ JS 施工引擎 = Python construction/launching/segmental（H1-H6）'));
process.exit(fail ? 1 : 0);
