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

  if (typeof module !== 'undefined' && module.exports) module.exports = BC;
  global.BC = BC;
})(typeof window !== 'undefined' ? window : globalThis);
