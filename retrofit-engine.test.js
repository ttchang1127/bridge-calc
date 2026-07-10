/* retrofit-engine.js 回歸測試：JS 引擎對 golden_answers.json 驗證（= Python bridgecalc.retrofit）。
 * 中國 JTG/T J22 單軌。執行：node retrofit-engine.test.js  （CI 亦跑；不符即非零退出）
 *
 * R1/R2/R4 同一原梁：b400 / h800 / h0=750 / As=1964(4D25) / C30（f_cd=13.8, f_sd=330）。
 * R2、R4 共用同一彈性開裂換算 x1（納入引擎時逼出 R2 算例 x1 修正：117→191.3）。
 */
var RF = require('./retrofit-engine.js');
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

console.log('\n' + pass + '/' + (pass + fail) + ' 對齊 golden' +
  (fail ? ' ❌ ' + fail + ' 項不符' : ' ✓ JS 補強引擎 = Python retrofit.py（R1/R2/R4）'));
process.exit(fail ? 1 : 0);
