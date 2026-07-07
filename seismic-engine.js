/* seismic-engine.js — seismic.py 的 JS 移植（台灣公路橋梁耐震設計規範，單軌）。
   對 golden seismic_S1/S2/S3/S5 驗證（雙實作不漂移）。
   沿用各條文原生單位（S1 長度 cm/L,H m；S2 力 kN/長度 m/勁度 kN/m；
   S3 kgf,cm²,tf·m；S5 無因次）。不做跨軌換算。 */
(function (global) {
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

  if (typeof module !== 'undefined' && module.exports) module.exports = SE;
  global.SE = SE;
})(typeof window !== 'undefined' ? window : globalThis);
