/* bridgecalc — JS 移植（與 Python engine/bridgecalc 同公式）
 * 簡支後張箱梁設計分析器核心。對 golden_answers.json 交叉驗證（台灣 HS20-44/8組×19股）。
 * 單位：力 N、長度 mm、應力 MPa；彎矩介面 kN·m（內部 ×1e6）；Pe 為 N。
 * node 可直接 require 做回歸；瀏覽器掛在 window.BC。
 */
(function (global) {
  'use strict';
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

  if (typeof module !== 'undefined' && module.exports) module.exports = BC;
  global.BC = BC;
})(typeof window !== 'undefined' ? window : globalThis);
