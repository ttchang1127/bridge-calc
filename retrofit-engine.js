/* retrofit-engine.js — retrofit.py 的 JS 移植（中國 JTG/T J22 補強承載力）。
   對 golden retrofit_* 驗證（雙實作不漂移）。單位 MPa/mm/N，彎矩回 kN·m。 */
(function (global) {
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

  if (typeof module !== 'undefined' && module.exports) module.exports = RF;
  global.RF = RF;
})(typeof window !== 'undefined' ? window : globalThis);
