/* seismic-engine.js 回歸測試：JS 引擎對 golden_answers.json 驗證（= Python bridgecalc.seismic）。
 * 台灣單軌（公路橋梁耐震設計規範）。執行：node seismic-engine.test.js  （CI 亦跑；不符即非零退出）
 */
var SE = require('./seismic-engine.js');
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

console.log('\n' + pass + '/' + (pass + fail) + ' 對齊 golden' +
  (fail ? ' ❌ ' + fail + ' 項不符' : ' ✓ JS 耐震引擎 = Python seismic.py（S1/S2/S3/S5）'));
process.exit(fail ? 1 : 0);
