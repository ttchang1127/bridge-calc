/* engine.js — 橋梁計算單一引擎（四域合一：箱梁 BC ＋ 耐震 SE ＋ 施工 CE ＋ 補強 RF）。
 * 由原四個 per-domain 引擎（box-girder/seismic/construction/retrofit-engine.js）收斂而成，
 * 消除「多檔各自與 Python 漂移」的面。瀏覽器掛 window.BC/SE/CE/RF（back-compat，呼叫端零改）；
 * node: const {BC,SE,CE,RF} = require('./engine.js')。對 golden 由 engine.test.js 一次驗 110 項。
 * 單一 closure → CE 直接用 BC.stresses（免原 global.BC 耦合）。
 */
(function (global) {
  'use strict';

  // ═══════════════ 箱梁主引擎（簡支箱梁核心＋次檢核＋連續梁）（原 box-girder-engine.js）═══════════════

  var BC = {};
  var sqrt = Math.sqrt, exp = Math.exp, abs = Math.abs, min = Math.min, max = Math.max;

  // ── 斷面 ──────────────────────────────────────────────
  BC.section = function (A, I, yb, h) {
    return { A: A, I: I, yb: yb, h: h, yt: h - yb, St: I / (h - yb), Sb: I / yb };
  };
  // 由構件尺寸算斷面性質（單室箱梁：頂板/底板/腹板）
  BC.sectionFromDims = function (topW, topT, botW, botT, webT, nWeb, h) {
    var hw = h - topT - botT;                 // 腹板淨高
    var p = [
      { A: topW * topT, y: h - topT / 2, I0: topW * Math.pow(topT, 3) / 12 },
      { A: botW * botT, y: botT / 2,     I0: botW * Math.pow(botT, 3) / 12 },
      { A: nWeb * webT * hw, y: botT + hw / 2, I0: nWeb * webT * Math.pow(hw, 3) / 12 }
    ];
    var A = 0, Ay = 0; p.forEach(function (s) { A += s.A; Ay += s.A * s.y; });
    var yb = Ay / A, I = 0;
    p.forEach(function (s) { I += s.I0 + s.A * Math.pow(s.y - yb, 2); });
    return BC.section(A, I, yb, h);
  };
  BC.tendon = function (n, strands, e, Ap, fpu, fpjRatio) {
    Ap = Ap || 140; fpu = fpu || 1860; fpjRatio = fpjRatio || 0.75;
    var Aps = n * strands * Ap, fpj = fpjRatio * fpu;
    return { n: n, strands: strands, e: e, Aps: Aps, fpu: fpu, fpj: fpj, Pi: fpj * Aps };
  };

  // ── 台灣 HS20-44 活載（卡車或車道取大）──────────────────
  BC.TW = { P: [36, 144, 144], x: [0, 4.25, 8.5], lane: 9.4, PM: 80, PV: 116 };
  function ilM(L, a, p) { return p <= a ? p * (L - a) / L : a * (L - p) / L; }
  BC.taiwanTruckMoment = function (L) {
    var T = BC.TW, best = 0;
    for (var a = 0.2; a < L; a += 0.2)
      for (var s = -T.x[2]; s <= L; s += 0.2) {
        var t = 0; for (var i = 0; i < 3; i++) { var xx = s + T.x[i]; if (xx >= 0 && xx <= L) t += T.P[i] * ilM(L, a, xx); }
        if (abs(t) > abs(best)) best = t;
      }
    return best;
  };
  BC.taiwanImpact = function (L) { return min(0.30, 15.24 / (L + 38.1)); };
  BC.taiwanLaneMoment = function (L) { return BC.TW.lane * L * L / 8 + BC.TW.PM * L / 4; };
  BC.taiwanPerLaneMoment = function (L) { return max(BC.taiwanTruckMoment(L), BC.taiwanLaneMoment(L)) * (1 + BC.taiwanImpact(L)); };
  function ilV(L, a, p) { return p < a ? -p / L : (L - p) / L; }
  BC.taiwanTruckShear = function (L) {     // 支承最大剪力（a→0）
    var T = BC.TW, best = 0, a = 1e-6;
    for (var s = -T.x[2]; s <= L; s += 0.1) {
      var t = 0; for (var i = 0; i < 3; i++) { var xx = s + T.x[i]; if (xx >= 0 && xx <= L) t += T.P[i] * ilV(L, a, xx); }
      if (abs(t) > abs(best)) best = t;
    }
    return abs(best);
  };
  BC.taiwanLaneShear = function (L) { return BC.TW.lane * L / 2 + BC.TW.PV; };
  BC.taiwanPerLaneShear = function (L) { return max(BC.taiwanTruckShear(L), BC.taiwanLaneShear(L)) * (1 + BC.taiwanImpact(L)); };
  BC.laneLiveLoad = function (perLaneM, nLanes, m) { return perLaneM * nLanes * (m == null ? 1 : m); };

  // ── 損失（摩擦/ES/潛變/乾縮/鬆弛）─────────────────────────
  BC.computeLosses = function (t, sec, M_DC, M_DW, o) {
    o = o || {};
    var mu = o.mu == null ? 0.25 : o.mu, K = o.K == null ? 0.003 : o.K,
        alpha = o.alpha == null ? 0.111 : o.alpha, x = o.x_ctrl == null ? 20 : o.x_ctrl,
        RH = o.RH == null ? 75 : o.RH, relax = o.relax == null ? 10 : o.relax,
        EpEci = o.Ep_Eci == null ? 7.33 : o.Ep_Eci;
    var Pi = t.Pi, e = t.e, fcgp = Pi / sec.A + Pi * e * e / sec.I - M_DC * 1e6 * e / sec.I;
    var fcds = M_DW * 1e6 * e / sec.I;
    var friction = t.fpj * (1 - exp(-(K * x + mu * alpha)));
    var shrink = 0.8 * (1195 - 10.55 * RH) * 0.0981;
    var ES = (t.n - 1) / (2 * t.n) * EpEci * fcgp;
    var creep = 12 * fcgp - 7 * fcds;
    var total = friction + ES + creep + shrink + relax;
    var fpe = t.fpj - total;
    return { fcgp: fcgp, ES: ES, creep: creep, friction: friction, shrink: shrink, relax: relax,
             total: total, loss_pct: total / t.fpj, fpe: fpe, Pe: fpe * t.Aps };
  };

  // ── 載重組合 ──────────────────────────────────────────
  BC.combinations = function (M_DC, M_DW, M_LL) {
    return { Strength_I: 1.25 * M_DC + 1.50 * M_DW + 1.75 * M_LL,
             Service_I: M_DC + M_DW + M_LL, Service_III: M_DC + M_DW + 0.8 * M_LL };
  };

  // ── 服務性應力 C1（壓為負）+ 最小預力反解 ─────────────────
  BC.stresses = function (Pe, sec, e, M_kNm) {
    var M = M_kNm * 1e6;
    return { st: -Pe / sec.A + Pe * e / sec.St - M / sec.St,
             sb: -Pe / sec.A - Pe * e / sec.Sb + M / sec.Sb };
  };
  BC.PeMinZeroTension = function (sec, e, M_kNm) { return M_kNm * 1e6 * sec.A / (sec.Sb + sec.A * e); };

  // ── 極限彎曲 M1 ───────────────────────────────────────
  BC.beta1 = function (fc) { return fc <= 28 ? 0.85 : max(0.65, 0.85 - 0.05 * (fc - 28) / 7); };
  BC.flexuralStrength = function (t, sec, fc, bEff, hf, dp, Mu, Pe, e, dt) {
    var Aps = t.Aps, fpu = t.fpu, fpy = 0.90 * fpu, k = 2 * (1.04 - fpy / fpu), b1 = BC.beta1(fc);
    var c = Aps * fpu / (0.85 * fc * b1 * bEff + k * Aps * fpu / dp);
    var a = b1 * c, fps = fpu * (1 - k * c / dp), Mn = Aps * fps * (dp - a / 2) / 1e6;
    dt = dt == null ? dp : dt;
    var eps_t = (dt - c) / c * 0.003, phi = eps_t >= 0.005 ? 1.0 : eps_t <= 0.002 ? 0.75 : 0.75 + 0.25 * (eps_t - 0.002) / 0.003;
    var phiMn = phi * Mn, fr = 0.97 * sqrt(fc), fcpe = Pe / sec.A + Pe * e / sec.Sb;
    var Mcr = sec.Sb * (fr + fcpe) / 1e6, lower = min(1.33 * Mu, 1.2 * Mcr);
    return { c: c, in_flange: c <= hf, a: a, fps: fps, Mn: Mn, eps_t: eps_t, phi: phi,
             phiMn: phiMn, CR: phiMn / Mu, Mcr: Mcr, lower: lower, ok: phiMn >= Mu && phiMn >= lower };
  };

  // ── 腹板抗剪 D1 ───────────────────────────────────────
  BC.principalTensionLimitTW = function (fc) { return 0.094 * sqrt(fc); };
  BC.shearWeb = function (Pe, sec, e, fc, bw, dv, Vu, xCtrl, L, nWebs, phi, fsy) {
    nWebs = nWebs || 2; phi = phi || 0.85; fsy = fsy || 420;
    var fpc = Pe / sec.A, slope = (4 * e / L) * (1 - 2 * xCtrl / L), Vp = (Pe / nWebs) * slope;
    var tau = Vu / (bw * dv), sigma1 = -fpc / 2 + sqrt(Math.pow(fpc / 2, 2) + tau * tau);
    var lim = BC.principalTensionLimitTW(fc), vc = 0.094 * sqrt(fc);
    var Vcw = (vc + 0.3 * fpc) * bw * dv + Vp, Vs_req = (Vu - phi * Vcw) / phi;
    return { fpc: fpc, slope: slope, Vp: Vp, tau: tau, sigma1: sigma1, sigma1_limit: lim,
             sigma1_ok: sigma1 <= lim, Vcw: Vcw, Vs_req: Vs_req, Av_s_req: max(Vs_req, 0) / (fsy * dv) };
  };
  BC.phiVn = function (Vcw, Av_s, dv, phi, fsy) { phi = phi || 0.85; fsy = fsy || 420; return phi * (Vcw + Av_s * fsy * dv); };
  BC.AvSminTW = function (fc, bw, fsy) { fsy = fsy || 420; return max(0.2 * sqrt(fc) * bw / fsy, 0.35 * bw / fsy); };

  // ── 撓度/預拱 C2/C3 ───────────────────────────────────
  BC.deflection = function (L, Ec, sec, w_DL, Pe, e, w_LL, phi, settle) {
    phi = phi == null ? 2.0 : phi; settle = settle == null ? 5.0 : settle;
    var K = 5 * Math.pow(L, 4) / (384 * Ec * sec.I), d_DL = K * w_DL, w_eq = 8 * Pe * e / (L * L);
    var d_PT = K * w_eq, net_el = d_DL - d_PT, net_LT = net_el * (1 + phi), d_LL = K * w_LL;
    return { K: K, d_DL: d_DL, w_eq: w_eq, d_PT: d_PT, net_elastic: net_el, LBR: d_PT / d_DL,
             net_long_term: net_LT, d_LL: d_LL, camber: max(net_LT, 0) + settle, d_LL_ok: d_LL <= L / 800 };
  };

  // ── 扭力 D2 ───────────────────────────────────────────
  BC.torsionCheck = function (Pe, sec, fc, Acp, pc, Tu_kNm, phi) {
    phi = phi || 0.90;
    var fpc = Pe / sec.A, k = 0.125 * sqrt(fc);
    var Tcr = k * (Acp * Acp / pc) * sqrt(1 + fpc / k) / 1e6, threshold = 0.25 * phi * Tcr;
    return { fpc: fpc, Tcr: Tcr, threshold: threshold, neglect: Tu_kNm < threshold, need_closed_stirrup: true };
  };

  // ── 橫向 D3（RC 板撓曲，每公尺寬）─────────────────────────
  BC.slabFlexure = function (Mu, As, d, fc, fy, b, phi) {
    b = b || 1000; phi = phi || 0.9;
    var a = As * fy / (0.85 * fc * b), phiMn = phi * As * fy * (d - a / 2) / 1e6;
    return { a: a, phiMn: phiMn, ok: phiMn >= Mu };
  };

  // ── 支承 E1 ───────────────────────────────────────────
  BC.bearingCheck = function (Rmax, Rmin, R_LL, delta, h_rt, Lb, Wb, te, Gkgf, gLim) {
    te = te || 10; Gkgf = Gkgf || 8; gLim = gLim || 0.50;
    var A = Lb * Wb, gamma = delta / h_rt, sigTL = Rmax * 1e3 / A;
    var S = Lb * Wb / (2 * te * (Lb + Wb)), limKgf = min(112, 1.66 * Gkgf * S), sigLim = limKgf * 0.0981;
    var Hm = Gkgf * 0.0981 * A * delta / h_rt / 1e3, stab = min(Lb, Wb) / 3;
    return { gamma_s: gamma, sigma_TL: sigTL, shape_S: S, sigma_TL_limit: sigLim, H_m: Hm,
             gamma_ok: gamma <= gLim, sigma_ok: sigTL <= sigLim, stability_ok: h_rt <= stab,
             H_ok: Hm <= Rmin / 5, no_uplift: (Rmin - R_LL) > 0 };
  };

  // ── 伸縮縫 E2 ─────────────────────────────────────────
  BC.expansionJoint = function (dl_T, dl_c, dl_s, g_install, margin) {
    margin = margin || 1.05;
    var sh = dl_T + dl_c + dl_s, gmax = g_install + sh, cap = Math.round(gmax * margin);
    return { shortening: sh, g_max: gmax, capacity: cap, joint: cap <= 75 ? 'Strip Seal 75mm' : cap <= 100 ? 'Modular ≤100mm' : 'Modular >100mm' };
  };

  // ── 錨碇 F1 ───────────────────────────────────────────
  BC.anchorageCheck = function (Pi_total, n, a_plate, h_diaph, n_per_web, fy, lf) {
    fy = fy || 400; lf = lf || 1.2;
    var Pu = lf * Pi_total / n, Tb = 0.25 * Pu * (1 - a_plate / h_diaph);
    var Fspall = 0.02 * n_per_web * Pu;
    return { Pu: Pu, Tburst: Tb, sum_Tburst: n_per_web * Tb, As_burst: 1.2 * Tb * 1e3 / (0.85 * fy),
             Fspall: Fspall, As_spall: Fspall * 1e3 / (0.6 * fy) };
  };
  BC.spiralLocalBearing = function (Pu, Pult0, flat, A_core, s, Dsp) {
    var Pult = Pult0 + 4.1 * flat * A_core * Math.pow(1 - s / Dsp, 2) / 1e3;
    return { Pult: Pult, margin: Pult / Pu, ok: Pult >= Pu };
  };

  // ── 疲勞 P1 ───────────────────────────────────────────
  BC.fatigueCheck = function (sec, Pe, e, M_perm, dM_fat, fc, EpEc, gamma) {
    EpEc = EpEc || 6.6; gamma = gamma || 1.75;
    var dM = gamma * dM_fat * 1e6, dsig_ps = EpEc * (dM / sec.I) * e;
    var sig_perm_top = -Pe / sec.A + Pe * e / sec.St - M_perm * 1e6 / sec.St;
    var sig_c_max = abs(sig_perm_top) + dM / sec.I * sec.yt, lim = 0.40 * fc;
    return { dsig_ps: dsig_ps, sig_c_max: sig_c_max, ps_ok: dsig_ps <= 125, c_ok: sig_c_max <= lim };
  };
  BC.stirrupFatigue = function (dV, s, Av, dv, lim) {
    lim = lim || 165; var d = dV * 1e3 * s / (Av * dv); return { dfsv: d, ok: d <= lim };
  };

  // ── 溫度梯度自平衡應力 T1 ─────────────────────────────────
  // bands: [[y_top,y_bot,area,T_mean],...]；fibers: [[name,y,T],...]
  BC.selfEquilibratingStress = function (bands, I, yt, h, fibers, Ec, alpha, neg) {
    Ec = Ec || 30590; alpha = alpha || 1.08e-5; neg = neg == null ? -0.30 : neg;
    var At = 0, sTA = 0, sM = 0;
    bands.forEach(function (b) { var yc = (b[0] + b[1]) / 2; At += b[2]; sTA += b[3] * b[2]; sM += b[3] * (yt - yc) * b[2]; });
    var Tu = sTA / At, TL = (h / I) * sM, Ea = Ec * alpha, pos = {}, ng = {};
    fibers.forEach(function (f) { var Tse = f[2] - (Tu + TL * (yt - f[1]) / h), sg = -Ea * Tse; pos[f[0]] = sg; ng[f[0]] = neg * sg; });
    return { Tu: Tu, TL: TL, sigma_pos: pos, sigma_neg: ng };
  };
  BC.thermalServiceCheck = function (sig_thermal, sig_base, gTG) {
    gTG = gTG || 0.5; var t = sig_base + gTG * sig_thermal; return { total: t, ok: t <= 0 };
  };
  // 由箱梁構件尺寸建 AASHTO 正梯度溫度分層（T1=18@頂、5@300mm、線性歸零@400）
  BC.thermalBandsFromDims = function (topW, topT, botW, botT, webT, nWeb, h) {
    var Tat = function (y) { return y <= 300 ? 18 - 13 * y / 300 : y <= 400 ? 5 * (400 - y) / 100 : 0; };
    var web = nWeb * webT, yb1 = h - botT;
    var bands = [[0, topT, topW * topT, (Tat(0) + Tat(topT)) / 2]];
    if (topT < 300) bands.push([topT, 300, web * (300 - topT), (Tat(topT) + 5) / 2]);
    bands.push([300, 400, web * 100, 2.5]);
    bands.push([400, yb1, web * (yb1 - 400), 0]);
    bands.push([yb1, h, botW * botT, 0]);
    return bands;
  };

  // ── 連續梁中墩（次彎矩 M2 + T 斷面極限彎曲）──────────────
  BC.secondaryMoment = function (M_total, M1) { return M_total - M1; };
  BC.primaryMoment = function (layers) { return layers.reduce(function (s, L) { return s + L[0] * L[1]; }, 0); };
  BC.flexuralStrengthT = function (Aps, fpu, fc, b, hf, bw, dp, Mu, dt) {
    var fpy = 0.90 * fpu, k = 2 * (1.04 - fpy / fpu), b1 = BC.beta1(fc);
    var c_rect = Aps * fpu / (0.85 * fc * b1 * b + k * Aps * fpu / dp), flanged, c, a, fps, Mn;
    if (c_rect <= hf) { flanged = false; c = c_rect; a = b1 * c; fps = fpu * (1 - k * c / dp); Mn = Aps * fps * (dp - a / 2) / 1e6; }
    else { flanged = true; c = (Aps * fpu - 0.85 * fc * (b - bw) * hf) / (0.85 * fc * b1 * bw + k * Aps * fpu / dp); a = b1 * c; fps = fpu * (1 - k * c / dp); Mn = (Aps * fps * (dp - a / 2) + 0.85 * fc * (b - bw) * hf * (a / 2 - hf / 2)) / 1e6; }
    dt = dt == null ? dp : dt;
    var eps_t = (dt - c) / c * 0.003, phi = eps_t >= 0.005 ? 1.0 : eps_t <= 0.002 ? 0.75 : 0.75 + 0.25 * (eps_t - 0.002) / 0.003;
    return { c: c, flanged: flanged, a: a, fps: fps, Mn: Mn, eps_t: eps_t, phi: phi, phiMn: phi * Mn, CR: phi * Mn / Mu, ok: phi * Mn >= Mu };
  };

  // ═══════════════ 耐震（台灣單軌 S1/S2/S3/S5）（原 seismic-engine.js）═══════════════

  var SE = {};

  // ── S1 落橋防止（§8.5） ──
  // 地盤變位係數 ε_G（表§8.5）
  SE.EPSILON_G = { "第一類": 0.0025, "第二類": 0.00375, "第三類": 0.005, "臺北盆地": 0.00625 };

  // 最小梁端防落長度 min N_L（式8-10）[cm]：(50+0.25L+H)(1+S²/8000)
  SE.minFalloffLength = function (L, H, S) {
    S = S || 0;
    return (50 + 0.25 * L + 1.0 * H) * (1 + S * S / 8000);
  };
  // 地盤水平相對變位 u_G [cm]：ε_G·L_e·(S_III/S_II)
  SE.groundRelativeDisplacement = function (soilClass, Le, ratioSIIISII) {
    return SE.EPSILON_G[soilClass] * Le * ratioSIIISII;
  };
  // 防落長度需求 L_N [cm]：活動 max(minNL, u_R+u_G)；固定 minNL
  SE.requiredFalloffLength = function (minNL, uR, uG, movable) {
    uR = uR || 0; uG = uG || 0;
    if (movable === false) return minNL;
    return Math.max(minNL, uR + uG);
  };
  // 落橋防止裝置設計降伏強度下限 F_y = 1.5·R_d [kN]
  SE.restrainerYieldStrength = function (Rd) { return 1.5 * Rd; };

  // ── S2 隔震與消能（第7章） ──
  // LRB 雙線性：割線勁度 K_eff [kN/m]、等效阻尼比 ξ_eq
  SE.isolationBilinear = function (Qd, Kd, Dd) {
    var Keff = Qd / Dd + Kd;
    var ED = 4 * Qd * Dd;
    var xiEq = ED / (2 * Math.PI * Keff * Dd * Dd);
    return { Keff: Keff, xiEq: xiEq };
  };
  // 一秒週期阻尼修正 B_1（表3-1 線性內插）。xi 為分數。
  SE.dampingCorrectionB1 = function (xi) {
    var p = xi * 100.0;
    if (p <= 2) return 0.80;
    if (p <= 5) return 0.80 + (1.00 - 0.80) * (p - 2) / (5 - 2);
    if (p <= 10) return 1.00 + (1.25 - 1.00) * (p - 5) / (10 - 5);
    if (p <= 20) return 1.25 + (1.50 - 1.25) * (p - 10) / (20 - 10);
    return 1.50;
  };
  // 短週期阻尼修正 B_S（表3-1 線性內插）
  SE.dampingCorrectionBS = function (xi) {
    var p = xi * 100.0;
    if (p <= 2) return 0.80;
    if (p <= 5) return 0.80 + (1.00 - 0.80) * (p - 2) / (5 - 2);
    if (p <= 10) return 1.00 + (1.33 - 1.00) * (p - 5) / (10 - 5);
    if (p <= 20) return 1.33 + (1.60 - 1.33) * (p - 10) / (20 - 10);
    return 1.60;
  };
  // 隔震橋等效線性化迭代（解說 C7.3）。回傳收斂結果。
  // W[kN],Qd[kN],Kd[kN/m],S_II_1；opts:{Kp(None=剛性墩),T0_II,S_II_S,taipeiBasin,g,tol,maxIter,sD0}
  SE.isolationDesign = function (W, Qd, Kd, S_II_1, opts) {
    opts = opts || {};
    var Kp = opts.Kp == null ? null : opts.Kp;
    var taipeiBasin = !!opts.taipeiBasin, T0_II = opts.T0_II, S_II_S = opts.S_II_S;
    var g = opts.g == null ? 9.81 : opts.g;
    var tol = opts.tol == null ? 1e-4 : opts.tol;
    var maxIter = opts.maxIter == null ? 100 : opts.maxIter;
    var sD = opts.sD0 == null ? 0.20 : opts.sD0;
    var it = 0, Keff = 0, xiE = 0, Te = 0, B1 = 0, Sa = 0, Dd = 0;
    for (it = 1; it <= maxIter; it++) {
      var Ke, bl;
      if (Kp === null) {
        Dd = sD;
        bl = SE.isolationBilinear(Qd, Kd, Dd);
        Keff = bl.Keff; Ke = Keff; xiE = bl.xiEq;
      } else {
        var DP = (Qd + Kd * sD) / (Kp + Kd);   // 墩頂位移（式C7-2）
        Dd = sD - DP;                           // 隔震器位移（式C7-1）
        bl = SE.isolationBilinear(Qd, Kd, Dd);
        Keff = bl.Keff; Ke = Keff * Kp / (Keff + Kp); xiE = bl.xiEq;  // 串聯（式7-3）
      }
      Te = 2 * Math.PI * Math.sqrt(W / (g * Ke));
      B1 = SE.dampingCorrectionB1(xiE);
      Sa = taipeiBasin ? (T0_II * S_II_S / (B1 * Te))   // 式7-1c
                       : (S_II_1 / (B1 * Te));          // 式7-1b
      var newSD = Sa * Te * Te / (4 * Math.PI * Math.PI) * g;   // 式7-1a
      if (Math.abs(newSD - sD) < tol) { sD = newSD; break; }
      sD = newSD;
    }
    Dd = Kp === null ? sD : sD - (Qd + Kd * sD) / (Kp + Kd);
    return {
      Dd: Dd, Te: Te, Keff: Keff, xiE: xiE, B1: B1, Sa: Sa,
      Vb_secant: Keff * Dd, Vb_bilinear: Qd + Kd * Dd, iterations: it > maxIter ? maxIter : it
    };
  };

  // ── S3 橋墩韌性（§4.2/5.3） ──
  // 超強彎矩 M_p = 1.3·M_n（§4.2.1）
  SE.overstrengthMoment = function (Mn, overstrength) {
    return (overstrength == null ? 1.3 : overstrength) * Mn;
  };
  // 容量設計剪力 V_u = ΣM_p / L_c
  SE.capacityShear = function (sumMp, Lc) { return sumMp / Lc; };
  // 圓柱塑鉸區螺箍體積比 ρ_s = max(式5-5, 式5-6)
  SE.rhoSCircular = function (fc, fyh, Ag, Ac, Pe) {
    var r55 = 0.45 * (fc / fyh) * (Ag / Ac - 1);
    var r56 = 0.12 * (fc / fyh) * (0.5 + 1.25 * Pe / (fc * Ag));
    return Math.max(r55, r56);
  };
  // 矩柱橫箍總斷面積 A_sh = max(式5-7, 式5-8)[cm²]
  SE.AshRectangular = function (a, hc, fc, fyh, Ag, Ac, Pe) {
    var a57 = 0.30 * a * hc * (fc / fyh) * (Ag / Ac - 1);
    var a58 = 0.12 * a * hc * (fc / fyh) * (0.5 + 1.25 * Pe / (fc * Ag));
    return Math.max(a57, a58);
  };
  // 圍束筋垂直間距上限 a ≤ min(15, 短邊/4, 6d_b)[cm]
  SE.confinementSpacingLimit = function (shortSide, db) {
    return Math.min(15.0, shortSide / 4.0, 6.0 * db);
  };
  // 塑鉸區配置長度 ℓ0 = max(柱深, ℓc/6, 45)[cm]
  SE.plasticHingeLength = function (colDepth, Lc, floor) {
    return Math.max(colDepth, Lc / 6.0, floor == null ? 45.0 : floor);
  };

  // ── S5 液狀化（§8.1，表8-1） ──
  // 土壤參數折減係數 D_E。F_L≥1 或 x>20m → 1.0
  SE.liquefactionReductionDE = function (FL, x, Rs) {
    if (FL >= 1.0 || x > 20) return 1.0;
    var shallow = x <= 10, dense = Rs > 0.3;
    if (FL <= 1.0 / 3.0) {                    // 第一級
      if (shallow) return dense ? 1.0 / 6.0 : 0.0;
      return 1.0 / 3.0;
    }
    if (FL <= 2.0 / 3.0) {                    // 第二級
      if (shallow) return dense ? 2.0 / 3.0 : 1.0 / 3.0;
      return 2.0 / 3.0;
    }
    if (shallow) return dense ? 1.0 : 2.0 / 3.0;   // 第三級
    return 1.0;
  };

  // ═══════════════ 施工（H1-H6 施拉/懸臂/推進/節塊）（原 construction-engine.js）═══════════════

  var CE = {};

  CE.JOINT_MIN_COMPRESSION_MPa = 0.21;   // 拼裝期接縫最小壓應力（SSOT）
  CE.BONDED_PT_MIN_RATIO = 0.30;         // 內置黏結 PT 最小比例

  // ── H1/H2 支架/托架施拉應力歷程 ──
  // 施拉容許拉 0.25√f'ci [MPa]
  CE.transferTensionLimit = function (fci) { return 0.25 * Math.sqrt(fci); };
  // 施拉容許壓（壓為負）：一般 0.55f'ci / 節塊 0.60f'ci
  CE.transferCompLimit = function (fci, bridgeType) {
    return -(bridgeType === '節塊' ? 0.60 : 0.55) * fci;
  };
  // 施工階段跨中頂/底緣應力與判定。P[N]、M_sw[kN·m]（支架上≈0，脫架=M_DC）
  CE.stageStress = function (P, sec, e, MswKNm, fci, bridgeType) {
    var s = BC.stresses(P, sec, e, MswKNm);
    var tl = CE.transferTensionLimit(fci), cl = CE.transferCompLimit(fci, bridgeType);
    return { st: s.st, sb: s.sb, top_ok: s.st <= tl, bot_ok: s.sb >= cl, tLimit: tl, cLimit: cl };
  };
  // 支架上分批張拉 n_batch/n_total 組，自重未活化（M_sw 預設 0）
  CE.batchedTransfer = function (Pi, nBatch, nTotal, sec, e, fci, MswKNm, bridgeType) {
    return CE.stageStress(Pi * nBatch / nTotal, sec, e, MswKNm || 0, fci, bridgeType);
  };

  // ── H3 平衡懸臂 ──
  // 變深箱梁拋物線斷面高 h(x) = h_mid + (h_pier − h_mid)(x/半跨)²
  CE.variableDepth = function (xFromMid, hPier, hMid, halfSpan) {
    return hMid + (hPier - hMid) * Math.pow(xFromMid / halfSpan, 2);
  };
  // 逐步懸臂彎矩 = Σ(G_i·arm_i) + 掛籃 G_FT·arm_FT
  CE.cantileverMoment = function (weights, arms, ftLoad, ftArm) {
    ftLoad = ftLoad || 0; ftArm = ftArm || 0;
    var s = 0;
    for (var i = 0; i < weights.length; i++) s += weights[i] * arms[i];
    return s + ftLoad * ftArm;
  };
  // 長期下撓 ≈ δ_elastic·(1+φ)
  CE.longTermDeflection = function (dElastic, phi) { return dElastic * (1 + phi); };

  // ── H4 推進 ILM ──
  // 懸臂根部最大負彎矩 M⁻ = −w·Lc²/2 [kN·m]
  CE.launchingCantileverMoment = function (wkNpm, LcM) { return -wkNpm * LcM * LcM / 2; };
  // 跨中最大正彎矩 M⁺ ≈ w·L²/8 [kN·m]
  CE.launchingSpanMoment = function (wkNpm, LM) { return wkNpm * LM * LM / 8; };
  // 臨時置中預力需求 Pc = (M⁺/Zb + σ_res)·A [kN]
  CE.centricPrestressRequired = function (MposKNm, ZbMm3, AMm2, sigmaRes) {
    sigmaRes = sigmaRes == null ? 1.5 : sigmaRes;
    return (MposKNm * 1e6 / ZbMm3 + sigmaRes) * AMm2 / 1e3;
  };
  // 最大正彎矩工況底緣應力 σ_bot = −Pc/A + M⁺/Zb [MPa，壓為負]
  CE.launchingBottomStress = function (PcKN, AMm2, MposKNm, ZbMm3) {
    return -PcKN * 1e3 / AMm2 + MposKNm * 1e6 / ZbMm3;
  };
  // 所需置中鋼腱束數 = ⌈Pc / 每束施拉力⌉
  CE.nTendons = function (PcKN, PperTendonKN) { return Math.ceil(PcKN / PperTendonKN); };
  // 頂推力 F = μ_s·W_total [kN]
  CE.jackingForce = function (muS, WtotalKN) { return muS * WtotalKN; };
  // 滑動支承局部支壓 σ_ba = R / A_bearing [MPa]
  CE.bearingStress = function (RkN, AbearingMm2) { return RkN * 1e3 / AbearingMm2; };

  // ── H5/H6 預鑄節塊 ──
  // 預鑄節塊自重 W = Ac·L_seg·γ [kN]
  CE.segmentWeight = function (AcM2, LsegM, gamma) {
    return AcM2 * LsegM * (gamma == null ? 25.0 : gamma);
  };
  // 接縫拼裝期最小臨時預力 = σ_min·Ac [kN]（Ac[mm²]）
  CE.jointMinPrestress = function (AcMm2, sigmaMin) {
    return (sigmaMin == null ? CE.JOINT_MIN_COMPRESSION_MPa : sigmaMin) * AcMm2 / 1e3;
  };
  // 接縫面均勻壓應力 σ = P/Ac [MPa]
  CE.jointCompression = function (PkN, AcMm2) { return PkN * 1e3 / AcMm2; };
  // 剪力鍵設計承載力 = V_fuk·ξ₁ξ₂Φ [kN/鍵]
  CE.shearKeyDesignCapacity = function (VfukKN, xiFactor) { return VfukKN * xiFactor; };
  // 剪力鍵驗核比 = (V_sd/n_keys)/V_key_design（≤1 通過）
  CE.shearKeyUtilization = function (VsdKN, nKeys, VkeyDesignKN) {
    return (VsdKN / nKeys) / VkeyDesignKN;
  };
  // 黏結 PT 比例 = P_bonded/P_total（≥ BONDED_PT_MIN_RATIO）
  CE.bondedPtRatio = function (PbondedKN, PtotalKN) { return PbondedKN / PtotalKN; };

  // ═══════════════ 補強（中國 JTG R1/R2/R4）（原 retrofit-engine.js）═══════════════

  var RF = {};

  // ── 共用：開裂換算斷面（二次受力） ──
  RF.crackedNaDepth = function (b, h0, As, alphaEs, AsComp, asComp) {
    AsComp = AsComp || 0; asComp = asComp || 0;
    var A1 = alphaEs * (As + AsComp) / b;
    var B1 = 2 * alphaEs * (As * h0 + AsComp * asComp) / b;
    return Math.sqrt(A1 * A1 + B1) - A1;
  };
  RF.crackedInertia = function (b, x1, As, h0, alphaEs, AsComp, asComp) {
    AsComp = AsComp || 0; asComp = asComp || 0;
    return b * Math.pow(x1, 3) / 3 + alphaEs * As * Math.pow(h0 - x1, 2) + alphaEs * AsComp * Math.pow(x1 - asComp, 2);
  };
  RF.initialConcreteStrain = function (Md1kNm, x1, Ec, Icr) {
    return Md1kNm * 1e6 * x1 / (Ec * Icr);
  };

  // ── R1 碳纖維 CFRP（式6-42/6-44） ──
  RF.cfrpKm1 = function (nf, Ef, tf) {
    var v = nf * Ef * tf;
    return v <= 214000 ? 1 - v / 428000 : 1070000 / v;
  };
  RF.cfrpAllowableStrain = function (nf, Ef, tf, epsFu, km2) {
    km2 = km2 == null ? 0.85 : km2;
    var km = Math.min(RF.cfrpKm1(nf, Ef, tf), km2, 0.9);
    return Math.min(km * epsFu, 2 / 3 * epsFu, 0.007);
  };
  RF.xiFb = function (epsCu, epsF, eps1) { return 0.8 * epsCu / (epsCu + epsF + eps1); };
  RF.cfrpMomentCapacity = function (b, h, h0, fcd, epsCu, As, fsd, Af, Ef, epsFAllow, eps1) {
    var xfb = RF.xiFb(epsCu, epsFAllow, eps1);
    var x = (fsd * As + Ef * epsFAllow * Af) / (fcd * b);
    var Mu = fsd * As * (h0 - 0.5 * xfb * h) + Ef * epsFAllow * Af * h * (1 - 0.5 * xfb);
    return { x: x, xi_fb: xfb, Mu_kNm: Mu / 1e6, case2: x <= xfb * h };
  };

  // ── R2 外貼鋼板（式6-26/6-35，ε_sp 迭代） ──
  RF.plateMomentCapacity = function (b, h, h0, fcd1, epsCu, As, fsd, Asp, fsp, Esp, x1, epsC1, beta, aS) {
    beta = beta || 0.8; if (aS == null) aS = h - h0;
    var sigmaSp = fsp, x = 0, epsSp = 0;
    for (var i = 0; i < 50; i++) {
      x = (fsd * As + sigmaSp * Asp) / (fcd1 * b);
      epsSp = epsCu * (beta * h - x) / x - epsC1 * (h - x1) / x1;
      var ns = Math.min(Esp * epsSp, fsp);
      if (Math.abs(ns - sigmaSp) < 0.01) { sigmaSp = ns; break; }
      sigmaSp = ns;
    }
    var Mu = fcd1 * b * x * (h0 - x / 2) + sigmaSp * Asp * aS;
    return { x: x, Mu_kNm: Mu / 1e6, plate_yields: Esp * epsSp >= fsp, sigma_sp: sigmaSp };
  };
  RF.plateDevLength = function (fsp, Asp, tauP, bP) { return fsp * Asp / (tauP * bP) + 300; };

  // ── R4 增大截面（式6-2/6-10，ε_s2 迭代） ──
  RF.enlargementMomentCapacity = function (b2, h02, h0, fcd1, epsCu, As1, fsd1, As2, fsd2, Es2, x1, epsC1, beta, AsComp, fsdComp, asComp) {
    beta = beta || 0.8; AsComp = AsComp || 0; fsdComp = fsdComp || 0; asComp = asComp || 0;
    var sigmaS2 = fsd2, x = 0, epsS2 = 0;
    for (var i = 0; i < 50; i++) {
      x = (fsd1 * As1 - fsdComp * AsComp + sigmaS2 * As2) / (fcd1 * b2);
      epsS2 = epsCu * (beta * h02 - x) / x - epsC1 * (h02 - x1) / x1;
      var ns = Math.min(Es2 * epsS2, fsd2);
      if (Math.abs(ns - sigmaS2) < 0.01) { sigmaS2 = ns; break; }
      sigmaS2 = ns;
    }
    var Mu = fcd1 * b2 * x * (h0 - x / 2) + fsdComp * AsComp * (h0 - asComp);
    return { x: x, Mu_kNm: Mu / 1e6, added_bar_yields: Es2 * epsS2 >= fsd2, sigma_s2: sigmaS2 };
  };

  // 原構件（單筋矩形）標稱抗彎 M0（供增幅對照）
  RF.baseMoment = function (b, h0, fcd, As, fsd) {
    var x = fsd * As / (fcd * b);
    return { x: x, M0_kNm: fsd * As * (h0 - x / 2) / 1e6 };
  };

  var API = { BC: BC, SE: SE, CE: CE, RF: RF };
  if (typeof module !== 'undefined' && module.exports) module.exports = API;
  global.BC = BC; global.SE = SE; global.CE = CE; global.RF = RF; global.Engine = API;
})(typeof window !== 'undefined' ? window : globalThis);
