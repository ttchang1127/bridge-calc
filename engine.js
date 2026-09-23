/* engine.js — 橋梁計算單一引擎（四域合一：箱梁 BC ＋ 耐震 SE ＋ 施工 CE ＋ 補強 RF）。
 * 由原四個 per-domain 引擎（box-girder/seismic/construction/retrofit-engine.js）收斂而成，
 * 消除「多檔各自與 Python 漂移」的面。瀏覽器掛 window.BC/SE/CE/RF（back-compat，呼叫端零改）；
 * node: const {BC,SE,CE,RF} = require('./engine.js')。對 golden 由 engine.test.js 一次驗 322 項。
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
  // 軸組雙向（原向＋掉頭）掃描取 |ΣP·η| 最大；整數步進、軸位夾到 [0,L]（同 influence.py _max_moving）
  function truckMaxMoving(L, il, step) {
    var T = BC.TW, span = T.x[2], n = Math.round((L + span) / step), best = 0;
    var sets = [[T.P, T.x], [T.P.slice().reverse(), T.x.map(function (d) { return span - d; }).reverse()]];
    for (var o = 0; o < 2; o++) for (var k = 0; k <= n; k++) {
      var s = -span + k * step, t = 0;
      for (var i = 0; i < 3; i++) { var p = s + sets[o][1][i]; if (p >= -1e-9 && p <= L + 1e-9) t += sets[o][0][i] * il(min(max(p, 0), L)); }
      if (abs(t) > abs(best)) best = t;
    }
    return best;
  }
  BC.taiwanTruckMoment = function (L) {
    var best = 0;
    for (var j = 1; j * 0.2 < L - 1e-9; j++) {
      var a = j * 0.2, m = truckMaxMoving(L, function (p) { return ilM(L, a, p); }, 0.2);
      if (abs(m) > abs(best)) best = m;
    }
    return best;
  };
  BC.taiwanImpact = function (L) { return min(0.30, 15.24 / (L + 38.1)); };
  BC.taiwanLaneMoment = function (L) { return BC.TW.lane * L * L / 8 + BC.TW.PM * L / 4; };
  BC.taiwanPerLaneMoment = function (L) { return max(BC.taiwanTruckMoment(L), BC.taiwanLaneMoment(L)) * (1 + BC.taiwanImpact(L)); };
  function ilV(L, a, p) { return p < a ? -p / L : (L - p) / L; }
  BC.taiwanTruckShear = function (L) {     // 支承最大剪力（a=0 即反力；L≥8.5 時 = 324−918/L）
    return abs(truckMaxMoving(L, function (p) { return ilV(L, 0, p); }, 0.1));
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
    // fricRatio 給定時直接採用（由 BC.frictionAt 依張拉端配置與線形算得），
    // 否則走原本的單點 α/x_ctrl 式——不傳即與既有結果完全相同。
    var fr = o.fricRatio == null ? (1 - exp(-(K * x + mu * alpha))) : o.fricRatio;
    var friction = t.fpj * fr;
    var shrink = 0.8 * (1195 - 10.55 * RH) * 0.0981;
    var ES = (t.n - 1) / (2 * t.n) * EpEci * fcgp;
    var creep = 12 * fcgp - 7 * fcds;
    var total = friction + ES + creep + shrink + relax;
    var fpe = t.fpj - total;
    return { fcgp: fcgp, ES: ES, creep: creep, friction: friction, shrink: shrink, relax: relax,
             total: total, loss_pct: total / t.fpj, fpe: fpe, Pe: fpe * t.Aps };
  };

  // ── 沿長度的預力 P(x)（同 prestress.loss_profile）──────────
  // computeLosses 是單點式：一個 fricRatio 代表全長。下面這支逐點算，
  // 且用**該點扣掉摩擦與滑移後的初始力**算 f_cgp（computeLosses 一律用 fpj·Aps）。
  // 差別最大的地方在錨碇端——摩擦≈0 但滑移最大，單點式兩者都取控制斷面值。
  BC.parabolicE = function (eMid, eEnd, L) {
    return function (x) { return eEnd + (eMid - eEnd) * 4 * x * (L - x) / (L * L); };
  };
  BC.udlMoment = function (w, L) { return function (x) { return w * x * (L - x) / 2; }; };

  BC.lossProfile = function (t, sec, L, eFn, MDFn, segs, o) {
    o = o || {};
    var mu = o.mu == null ? 0.25 : o.mu, K = o.K == null ? 0.003 : o.K,
        jack = o.jack || 'both', slip = o.slip == null ? 0 : o.slip,
        Ep = o.Ep == null ? 195000 : o.Ep, anchors = o.anchors || null,
        RH = o.RH == null ? 75 : o.RH, relax = o.relax == null ? 10 : o.relax,
        EpEci = o.Ep_Eci == null ? 7.33 : o.Ep_Eci, n = o.n == null ? 20 : o.n,
        MSDLFn = o.M_SDL_fn || null;
    var fpj = t.fpj, Aps = t.Aps, N = t.n;
    var shrink = 0.8 * (1195 - 10.55 * RH) * 0.0981;
    var xs = [], pts = [], i;
    for (i = 0; i <= n; i++) {
      var x = L * i / n, fric, sl;
      if (anchors) {
        var f = BC.segmentedTendonForce(x, anchors, segs, fpj, Aps, 0, mu, K, jack, slip, Ep);
        fric = f.friction; sl = f.slip;
      } else {
        fric = fpj * BC.frictionAt(x, L, segs, mu, K, jack);
        sl = BC.tendonSlipLoss(x, L, segs, fpj, jack, mu, K, slip, Ep);
      }
      var Pix = (fpj - fric - sl) * Aps, e = eFn(x),
          M_D = MDFn(x) * 1e6, M_SDL = MSDLFn ? MSDLFn(x) * 1e6 : 0,
          fcgp = Pix / sec.A + Pix * e * e / sec.I - M_D * e / sec.I,
          fcds = M_SDL * e / sec.I,
          ES = (N - 1) / (2 * N) * EpEci * fcgp,
          creep = 12 * fcgp - 7 * fcds,
          shortL = fric + sl + ES, lng = creep + shrink + relax,
          total = shortL + lng, fpe = fpj - total;
      xs.push(x);
      pts.push({ x: x, e: e, fric: fric, slip: sl, ES: ES, creep: creep, shrink: shrink,
                 relax: relax, short: shortL, long: lng, total: total, loss_pct: total / fpj,
                 fcgp: fcgp, fpe: fpe, Pe: fpe * Aps });
    }
    var Pes = pts.map(function (p) { return p.Pe; }),
        lps = pts.map(function (p) { return p.loss_pct; }),
        imin = 0, imax = 0, ilmax = 0, lmin = lps[0], area = 0;
    // 嚴格不等號 → 平手取最先者，與 Python list.index(min(...)) 一致（兩端對稱時取 x=0）
    for (i = 1; i <= n; i++) {
      if (Pes[i] < Pes[imin]) imin = i;
      if (Pes[i] > Pes[imax]) imax = i;
      if (lps[i] > lps[ilmax]) ilmax = i;
      if (lps[i] < lmin) lmin = lps[i];
    }
    for (i = 0; i < n; i++) area += (Pes[i] + Pes[i + 1]) / 2 * (xs[i + 1] - xs[i]);
    return {
      xs: xs, pts: pts, Pe_min: Pes[imin], x_Pemin: xs[imin], Pe_max: Pes[imax],
      x_Pemax: xs[imax], Pe_avg: L ? area / L : 0, loss_pct_max: lps[ilmax],
      x_lossmax: xs[ilmax], loss_pct_min: lmin,
      // 取樣點之間線性內插；f_cgp 沿 x 為二次，故非取樣點有內插誤差（調大 n 可減小）
      at: function (xq) {
        if (xq <= xs[0]) return pts[0];
        if (xq >= xs[xs.length - 1]) return pts[pts.length - 1];
        for (var j = 0; j < xs.length - 1; j++) {
          if (xs[j] <= xq && xq <= xs[j + 1]) {
            var tt = xs[j + 1] > xs[j] ? (xq - xs[j]) / (xs[j + 1] - xs[j]) : 0,
                a = pts[j], b = pts[j + 1],
                lp = function (u, v) { return u + (v - u) * tt; };
            return { x: xq, e: lp(a.e, b.e), fric: lp(a.fric, b.fric), slip: lp(a.slip, b.slip),
                     ES: lp(a.ES, b.ES), creep: lp(a.creep, b.creep), shrink: a.shrink,
                     relax: a.relax, short: lp(a.short, b.short), long: lp(a.long, b.long),
                     total: lp(a.total, b.total), loss_pct: lp(a.loss_pct, b.loss_pct),
                     fcgp: lp(a.fcgp, b.fcgp), fpe: lp(a.fpe, b.fpe), Pe: lp(a.Pe, b.Pe) };
          }
        }
        return pts[pts.length - 1];
      }
    };
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

  // ── 反解設計（design.py 移植）：設計＝驗算之逆 ───────────
  // Pe = n·股數·A_strand·f_pe → n = ⌈Pe/(股數·A_strand·f_pe)⌉
  BC.minTendonGroups = function (Pe_req_N, strandsPer, fpe, Ap) {
    Ap = Ap == null ? 140 : Ap;
    return Math.ceil(Pe_req_N / (strandsPer * Ap * fpe));
  };
  // w_eq = 8·Pe·a/L²、LBR = w_eq/w_DL → a = LBR·w_DL·L²/(8·Pe)（w_DL N/mm、L mm → a mm）
  BC.requiredDrape = function (LBR_target, w_DL, L, Pe) {
    return LBR_target * w_DL * L * L / (8 * Pe);
  };
  // σ_b = −Pe/A − Pe·e/S_b + M/S_b ≤ σ_lim → S_b ≥ (M − Pe·e)/(σ_lim + Pe/A)
  BC.minSectionModulusSb = function (Pe, A, e, M_kNm, sigma_tension_limit) {
    var lim = sigma_tension_limit == null ? 0 : sigma_tension_limit;
    return (M_kNm * 1e6 - Pe * e) / (lim + Pe / A);
  };

  // ── 摩擦損失沿長度分佈 G1（tendon_profile.friction_* 移植）──
  // ΔP/P = 1 − e^−(μ·α + K·x)：α 是自張拉端累積到該斷面的角變化、x 是該段路徑長，
  // 兩者都從張拉端起算 → 張拉端配置直接決定每個斷面的損失。
  // α(x)=∫|e″|dx，反曲時 e″ 變號必須逐段取絕對值累加（用斜率差會相消）。
  // 'both' 雙端：該點取損失小者；'alt' 交錯：半數腱自左半數自右 → 取兩者平均。
  // ⚠ 跨中在四種配置下同值（左右對稱、路徑各半），差異在其餘斷面。
  BC.parabolicCurvSegs = function (L, a) {     // 單跨拋物線 κ=|e″|=8a/(L²·1000) rad/m
    return [[0, L, 8 * a / (L * L * 1000)]];
  };
  BC.frictionAngle = function (segs, x, L, fromEnd) {
    var t = 0;
    for (var i = 0; i < segs.length; i++) {
      var a = segs[i][0], b = segs[i][1], k = segs[i][2];
      var lo = fromEnd ? Math.max(a, x) : Math.max(a, 0), hi = fromEnd ? Math.min(b, L) : Math.min(b, x);
      if (hi > lo) t += k * (hi - lo);
    }
    return t;
  };
  BC.frictionAt = function (x, L, segs, mu, K, jack) {
    mu = mu == null ? 0.25 : mu; K = K == null ? 0.003 : K; jack = jack || 'both';
    var lS = 1 - exp(-(mu * BC.frictionAngle(segs, x, L, false) + K * x));
    var lE = 1 - exp(-(mu * BC.frictionAngle(segs, x, L, true) + K * (L - x)));
    if (jack === 'start') return lS;
    if (jack === 'end') return lE;
    if (jack === 'both') return Math.min(lS, lE);
    return 0.5 * (lS + lE);
  };
  BC.frictionProfile = function (L, segs, mu, K, jack, n) {
    n = n || 20; jack = jack || 'both';
    var xs = [], rs = [], i, area = 0, imax = 0;
    for (i = 0; i <= n; i++) { xs.push(L * i / n); rs.push(BC.frictionAt(xs[i], L, segs, mu, K, jack)); }
    for (i = 0; i < n; i++) area += (rs[i] + rs[i + 1]) / 2 * (xs[i + 1] - xs[i]);
    for (i = 0; i <= n; i++) if (rs[i] > rs[imax]) imax = i;
    return { xs: xs, ratios: rs, at_mid: BC.frictionAt(L / 2, L, segs, mu, K, jack),
             at_start: rs[0], at_end: rs[n], avg: area / L,
             max_ratio: rs[imax], x_max: xs[imax], jack: jack };
  };

  // ── 逐腱合成 G1（tendon_profile.tendon_forces 移植）─────────
  // 總 Pe 用「平均損失率×總面積」其實等價（各項對 f_pe 線性）；真正的差別在**合力位置**：
  //   e_eff = Σ(Pe_i·e_i)/ΣPe_i ≠ e_geom（各腱 Pe 不同且高度不同時）
  // 交錯端張拉若落在不同層（單列多層＝序號就是層號），低層自起點、損失小、力大
  // → 合力下移；同層左右配對（多列）垂直相消，只剩橫向偏心 x_eff。
  BC.assignJack = function (n, mode) {
    var out = [], i;
    for (i = 0; i < n; i++) out.push(
      (mode === 'both' || mode === 'start' || mode === 'end') ? mode : (i % 2 === 0 ? 'start' : 'end'));
    return out;
  };
  // 單腱錨具滑移損失 MPa（同 tendon_profile.tendon_slip_loss）：'both' 兩端各影響半長、'start'/'end' 全長
  BC.tendonSlipLoss = function (x, L, segs, fpj, jack, mu, K, slip, Ep) {
    if (!slip || slip <= 0 || L <= 0) return 0;
    jack = jack || 'both';
    var Lm, d, r;
    // 雙端張拉：依 x 靠近哪一端取該端的摩擦梯度（線形不對稱時兩端不同）
    if (jack === 'both') { Lm = L / 2; if (x <= Lm) { d = x; r = BC.frictionAt(Lm, L, segs, mu, K, 'start'); } else { d = L - x; r = BC.frictionAt(Lm, L, segs, mu, K, 'end'); } }
    else if (jack === 'end') { Lm = L; d = L - x; r = BC.frictionAt(0, L, segs, mu, K, 'end'); }
    else { Lm = L; d = x; r = BC.frictionAt(L, L, segs, mu, K, 'start'); }
    return BC.anchorSlipLoss(d, Lm, Lm > 0 ? fpj * r / Lm : 0, slip, Ep).dsigma;
  };
  // 分段（中間錨碇）腱：摩擦與滑移只在該段內累積（同 tendon_profile.segmented_tendon_force）
  BC.segmentedTendonForce = function (x, anchors, curvSegs, fpj, Aps, otherLoss, mu, K, jack, slip, Ep) {
    jack = jack || 'both';
    for (var i = 0; i < anchors.length - 1; i++) {
      var x0 = anchors[i], x1 = anchors[i + 1];
      if (x >= x0 - 1e-9 && x <= x1 + 1e-9) {
        var L = x1 - x0, loc = [];
        curvSegs.forEach(function (c) { var lo = Math.max(c[0], x0), hi = Math.min(c[1], x1); if (hi > lo) loc.push([lo - x0, hi - x0, c[2]]); });
        var sx = Math.min(Math.max(x - x0, 0), L), r = BC.frictionAt(sx, L, loc, mu, K, jack);
        var sl = BC.tendonSlipLoss(sx, L, loc, fpj, jack, mu, K, slip, Ep), fpe = fpj - fpj * r - sl - otherLoss;
        return { P: fpe * Aps / 1000, fpe: fpe, friction: fpj * r, slip: sl, seg: i, L_seg: L, s: sx };
      }
    }
    return { P: 0, fpe: 0, friction: 0, slip: 0, seg: -1, L_seg: 0, s: 0 };
  };
  BC.segmentedFrictionProfile = function (anchors, curvSegs, fpj, mu, K, jack, slip, n) {
    n = n || 40;
    var tot = anchors[anchors.length - 1] - anchors[0], xs = [], rs = [], i;
    for (i = 0; i <= n; i++) {
      var x = anchors[0] + tot * i / n, f = BC.segmentedTendonForce(x, anchors, curvSegs, fpj, 1, 0, mu, K, jack, slip);
      xs.push(x); rs.push((f.friction + f.slip) / fpj);
    }
    var mx = Math.max.apply(null, rs);
    return { xs: xs, ratios: rs, max_ratio: mx, x_max: xs[rs.indexOf(mx)], avg: rs.reduce(function (a, b) { return a + b; }, 0) / rs.length };
  };
  BC.tendonForces = function (tendons, x, L, segs, yb, fpj, ApEach, otherLoss, mu, K, slip, Ep, anchors) {
    var per = [], sP = 0, sPe = 0, sPx = 0, sr = 0, i, sy = 0;
    for (i = 0; i < tendons.length; i++) {
      var t = tendons[i], jk = t.jack || 'both', r, sl, fpe;
      if (anchors && anchors.length > 2) { var sf = BC.segmentedTendonForce(x, anchors, segs, fpj, ApEach, otherLoss, mu, K, jk, slip, Ep);
        r = sf.friction / fpj; sl = sf.slip; fpe = sf.fpe; }
      else { r = BC.frictionAt(x, L, segs, mu, K, jk); sl = BC.tendonSlipLoss(x, L, segs, fpj, jk, mu, K, slip, Ep); fpe = fpj - fpj * r - sl - otherLoss; }
      var Pe = fpe * ApEach, tx = t.x || 0;
      per.push({ no: t.no || '', y: t.y, x: tx, jack: jk, ratio: r, slip: sl, fpe: fpe, Pe: Pe });
      sP += Pe; sPe += Pe * (yb - t.y); sPx += Pe * tx; sr += r; sy += t.y;
    }
    var n = tendons.length, eG = yb - sy / n, eE = sP ? sPe / sP : 0;
    return { Pe_total: sP, e_eff: eE, e_geom: eG, de: eE - eG, x_eff: sP ? sPx / sP : 0,
             Pe_avg_ratio: (fpj - fpj * (sr / n) - otherLoss) * ApEach * n, per: per };
  };

  // ── 管道實配排列 G1（tendon_profile.duct_layout 移植）───
  // 分析用的 CGS 是「n 組腱視為一點」的簡化；實配要排 n_col 列 × n_row 層，最底層
  // 必低於 CGS。排完整體平移使實配形心 = y_cgs，故引擎在用的偏心 e 不因實配改變。
  // 淨間距需求，兩種規則（⚠ 這個選擇直接決定腹板內排得下幾列）：
  //   'tw'（預設）max(40, 1.5·d_agg)——台灣 §8.25.2 明文兩條件，同公式卡 G1 §7.1 表列
  //   'od'         再取 max(…, 孔道外徑)——AASHTO「孔道 OD > 100 mm 時淨距 ≥ 孔道外徑」；
  //                2026-09-18 核對第八章原文 p.178–179 確認**非台灣規範要求**，保留為保守設計慣例。
  // φ100、腹板 350、保護層 40 時：'tw' → 兩列（需求 40）、'od' → 單列（需求 100）。
  BC.ductSpacingRequired = function (od, dAgg, rule) {
    dAgg = dAgg == null ? 25 : dAgg;
    var b = Math.max(40, 1.5 * dAgg);
    return rule === 'od' ? Math.max(b, od) : b;
  };
  // 拋物線最小曲率半徑 R = L²/(8a)（mm）
  // ⚠ 最小值 6,000（錨碇區 3,600、PE 管 9,000）是 **AASHTO 5.4.6.1**；日本 ≥100×管道外徑。
  //   台灣規範未規定最小曲率半徑（2026-09-20 NLM 核第八章，原註解誤標為台灣）。
  BC.radiusOfCurvature = function (a, L) { return L * L / (8 * a); };
  BC.ductLayout = function (nTendons, nWeb, yCgs, o) {
    o = o || {};
    var od = o.od == null ? 100 : o.od, webT = o.webT == null ? 350 : o.webT,
        cover = o.cover == null ? 40 : o.cover, sv = o.sv == null ? 40 : o.sv,
        dAgg = o.dAgg == null ? 25 : o.dAgg, yb = o.yb, h = o.h, rule = o.rule || 'tw',
        sideMargin = o.sideMargin || 0;   // 側向額外餘裕：預設 0 時最外側管保護層恰等於 cover
    nWeb = Math.max(1, nWeb | 0); var n = Math.max(1, nTendons | 0);
    var base = Math.floor(n / nWeb), rem = n % nWeb, perWeb = [], i, j;
    for (i = 0; i < nWeb; i++) perWeb.push(base + (i < rem ? 1 : 0));
    var nPerWeb = Math.max.apply(null, perWeb);

    var sReq = BC.ductSpacingRequired(od, dAgg, rule), avail = webT - 2 * (cover + sideMargin);
    var nCol = avail >= od ? Math.floor((avail + sReq) / (od + sReq)) : 0;
    nCol = Math.max(1, nCol);                       // 至少畫一列；放不下時由 sHOk 報 ✗
    var sH = nCol > 1 ? (avail - nCol * od) / (nCol - 1) : avail - od;
    var sHOk = avail >= od && (nCol === 1 || sH >= sReq - 1e-9);

    var nRow = Math.ceil(nPerWeb / nCol), pitchV = od + sv, counts = [];
    for (i = 0; i < nRow; i++) counts.push(Math.min(nCol, nPerWeb - i * nCol));
    var sum = 0; for (i = 0; i < nRow; i++) sum += counts[i] * i * pitchV;
    var relCgs = sum / nPerWeb, shift = yCgs - relCgs, rows = [];
    for (i = 0; i < nRow; i++) rows.push({ y: i * pitchV + shift, count: counts[i] });

    var pitchH = od + (nCol > 1 ? sH : 0), xs = [];
    for (j = 0; j < nCol; j++) xs.push((j - (nCol - 1) / 2) * pitchH);

    var yBot = rows[0].y, yTop = rows[nRow - 1].y, coverBot = yBot - od / 2;
    var coverOk = coverBot >= cover - 1e-9, svOk = sv >= sReq - 1e-9;
    var topOk = h == null ? true : (yTop + od / 2 <= h - cover + 1e-9);
    return { nPerWeb: nPerWeb, perWeb: perWeb, nCol: nCol, nRow: nRow, rows: rows,
             xOffsets: xs, pitchV: pitchV, sReq: sReq, sH: sH, sv: sv, webAvail: avail,
             yCgs: yCgs, yBot: yBot, yTop: yTop, coverBot: coverBot,
             coverOk: coverOk, svOk: svOk, sHOk: sHOk, topOk: topOk,
             fits: coverOk && svOk && sHOk && topOk,
             eMax: yb == null ? null : yb - (cover + od / 2 + relCgs),
             eMaxPoint: yb == null ? null : yb - cover - od / 2,
             coverSide: (webT - (xs[nCol - 1] - xs[0]) - od) / 2 };
  };

  // §8.25.3 套管捆紮（同 tendon_profile.duct_layout_bundled）：每束 ≤3 管接觸、束間 §8.25.2 淨距；
  // "H" 同層並排接觸（窄腹板決定性）、"V" 上下疊放接觸（壓低合力）。端部 90 cm 內須維持間距。
  BC.ductLayoutBundled = function (nTendons, nWeb, yCgs, o) {
    o = o || {};
    var base = BC.ductLayout(nTendons, nWeb, yCgs, o), od = o.od == null ? 100 : o.od,
        webT = o.webT == null ? 350 : o.webT, cover = o.cover == null ? 40 : o.cover,
        sv = o.sv == null ? 40 : o.sv, yb = o.yb, h = o.h, mb = Math.max(1, Math.min(3, o.maxBundle || 3));
    base.bundle = 'none';
    var cands = { none: base }, n = base.nPerWeb, avail = base.webAvail, i;
    function build(pos, tag, sH) {
      var rel = pos.reduce(function (a, p) { return a + p[1]; }, 0) / pos.length, shift = yCgs - rel;
      var ys = [], xs = [];
      pos.forEach(function (p) {
        if (!ys.some(function (v) { return Math.abs(v - p[1]) < 1e-9; })) ys.push(p[1]);
        if (!xs.some(function (v) { return Math.abs(v - p[0]) < 1e-9; })) xs.push(p[0]);
      });
      ys.sort(function (a, b) { return a - b; }); xs.sort(function (a, b) { return a - b; });
      var rows = ys.map(function (y) { return { y: y + shift, count: pos.filter(function (p) { return Math.abs(p[1] - y) < 1e-9; }).length }; });
      var yBot = rows[0].y, yTop = rows[rows.length - 1].y, coverBot = yBot - od / 2,
          coverOk = coverBot >= cover - 1e-9, topOk = h == null ? true : (yTop + od / 2 <= h - cover + 1e-9),
          svOk = sv >= base.sReq - 1e-9;
      return { nPerWeb: n, perWeb: base.perWeb, nCol: xs.length, nRow: rows.length, rows: rows, xOffsets: xs,
               pitchV: od + sv, sReq: base.sReq, sH: sH, sv: sv, webAvail: avail, yCgs: yCgs, yBot: yBot, yTop: yTop,
               coverBot: coverBot, coverOk: coverOk, svOk: svOk, sHOk: true, topOk: topOk, fits: coverOk && topOk && svOk,
               eMax: yb == null ? null : yb - (cover + od / 2 + rel), eMaxPoint: yb == null ? null : yb - cover - od / 2,
               coverSide: (webT - (xs[xs.length - 1] - xs[0]) - od) / 2, bundle: tag };
    }
    var k = avail >= od ? Math.min(mb, Math.floor(avail / od)) : 0;
    if (k > base.nCol) {
      var pH = [];
      for (i = 0; i < n; i++) { var r = Math.floor(i / k), j = i % k;   // 末層不另置中（與一般排列一致）
        pH.push([(j - (k - 1) / 2) * od, r * (od + sv)]); }
      cands.H = build(pH, 'H', 0);
    }
    if (base.nRow > 1 && mb > 1) {
      var pV = [], nc = base.nCol;
      for (i = 0; i < n; i++) { var layer = Math.floor(i / nc), c = i % nc, gg = Math.floor(layer / mb), jj = layer % mb;
        pV.push([c < base.xOffsets.length ? base.xOffsets[c] : 0, gg * (mb * od + sv) + jj * od]); }
      cands.V = build(pV, 'V', base.sH);
    }
    var best = base;
    Object.keys(cands).forEach(function (t) { var v = cands[t];
      if (v.fits && v.eMax != null && (!best.fits || v.eMax > best.eMax)) best = v; });
    return { best: best, cands: cands };
  };

  // 各管中心 [[x, y]]，依填充順序（＝編號 W<腹板>-<序號>）；一般/H/V 皆為第 i 管在 i/nCol 層、i%nCol 列
  BC.ductPositions = function (r) {
    var out = [];
    for (var i = 0; i < r.nPerWeb; i++) out.push([r.xOffsets[i % r.nCol], r.rows[Math.floor(i / r.nCol)].y]);
    return out;
  };
  // §8.25.3 端部 90 cm（同 tendon_profile.end_zone_duct_check）：距構材端部 90 cm 內須維持 §8.25.2 間距。
  // eValues＝端部範圍內取樣偏心；排列可行性對合力高度單調，只檢最小／最大偏心。bundle 非 'none' 時
  // 以端部邊界 CGS 比較捆紮與展開排列，逐管求過渡段移位。⚠ 移位造成的偏折力另檢。
  BC.endZoneDuctCheck = function (nTendons, nWeb, eValues, yb, h, o) {
    o = o || {};
    var bundle = o.bundle || 'none', endLen = o.endLen == null ? 900 : o.endLen,
        wt = o.webTEnd == null ? (o.webT == null ? 350 : o.webT) : o.webTEnd;
    var oe = Object.assign({}, o, { webT: wt, yb: yb, h: h }), ob = Object.assign({}, o, { yb: yb, h: h });
    var eLo = Math.min.apply(null, eValues), eHi = Math.max.apply(null, eValues);
    var lo = BC.ductLayout(nTendons, nWeb, yb - eLo, oe), hi = BC.ductLayout(nTendons, nWeb, yb - eHi, oe);
    var moves = [], dxMax = 0, dyMax = 0;
    if (bundle === 'H' || bundle === 'V') {
      var eB = eValues[eValues.length - 1], cands = BC.ductLayoutBundled(nTendons, nWeb, yb - eB, ob).cands;
      if (cands[bundle]) {
        var pb = BC.ductPositions(cands[bundle]), pu = BC.ductPositions(BC.ductLayout(nTendons, nWeb, yb - eB, oe));
        for (var i = 0; i < pb.length; i++) {
          var dx = pu[i][0] - pb[i][0], dy = pu[i][1] - pb[i][1];
          moves.push([i + 1, dx, dy]); dxMax = Math.max(dxMax, Math.abs(dx)); dyMax = Math.max(dyMax, Math.abs(dy));
        }
      }
    }
    return { endLen: endLen, webTEnd: wt, eLo: eLo, eHi: eHi, layLo: lo, layHi: hi, fits: lo.fits && hi.fits,
             bundle: bundle, dxMax: dxMax, dyMax: dyMax, moves: moves };
  };

  // 腹板厚與管徑（同 tendon_profile.web_duct_check）：台灣規範三項均未規定——§8.9.3 僅漸變長度 ≥ 12 倍厚度差；
  // 其餘為 AASHTO LRFD（2008 SI）參考：5.4.6.2 管徑 ≤ 0.4 × 該處最小厚度、C5.14.1.5.1c 最小腹板
  // 200/300/380、5.8.2.9 b_v 扣同高程管徑 ¼（灌漿）／½（未灌漿）。H 捆紮束寬比僅列出不判定。
  BC.AASHTO_WEB_MIN = { none: 200, one_way: 300, both: 380 };
  BC.webDuctCheck = function (webT, od, nLevel, h, verticalDucts, bundle, dtChange) {
    var ratio = od / webT, bw = bundle === 'H' ? nLevel * od : null,
        wmin = BC.AASHTO_WEB_MIN[verticalDucts ? 'both' : 'one_way'], s = nLevel * od;
    return { webT: webT, od: od, nLevel: nLevel, ratio: ratio, ratioOk: ratio <= 0.4 + 1e-12,
             bundleW: bw, bundleRatio: bw == null ? null : bw / webT, webMin: wmin, webMinOk: webT >= wmin - 1e-9,
             deepNote: h != null && h > 2400, bvGrouted: webT - 0.25 * s, bvUngrouted: webT - 0.5 * s,
             taperMin: dtChange == null ? null : 12 * Math.abs(dtChange) };
  };

  // 曲線鋼腱偏折力（同 tendon_profile.deviation_force_check）：AASHTO 5.10.4.3。
  // F_u−in / F_u−out = P_u/R（N/mm）；保護層抗拉脫 V_r = φ·0.33√f'ci·d_c、d_c = 保護層 + od/2；
  // 不足 → 面內配完全錨定 tie-back、面外沿曲線段配局部圍束（宜螺旋筋）。圍束 f_s ≤ 0.6f_y、f_y ≤ 420、
  // 間距 ≤ min(3·od, 600)。兩平面同時彎曲時向量相加。R ≥ 6,000（錨碇區 3,600、PE 9,000）＝5.4.6.1。
  // ⚠ 台灣規範未規定偏折力（§8.21.3 6.(2) 僅要求隔梁錨碇處配筋抵抗曲率分力）。
  BC.transitionRadius = function (offset, len) { return Math.abs(offset) < 1e-12 ? Infinity : len * len / (4 * Math.abs(offset)); };
  BC.deviationForceCheck = function (Pu, dx, dy, Lt, fci, coverSide, coverFace, od, phi) {
    od = od == null ? 100 : od; phi = phi == null ? 0.9 : phi;
    var Rl = BC.transitionRadius(dx, Lt), Rv = BC.transitionRadius(dy, Lt);
    var Fo = Rl === Infinity ? 0 : Pu / Rl, Fi = Rv === Infinity ? 0 : Pu / Rv;
    var dcl = coverSide + od / 2, dcv = coverFace + od / 2,
        vn = function (dc) { return phi * 0.33 * Math.sqrt(fci) * dc; }, Vl = vn(dcl), Vv = vn(dcv);
    var ltm = function (d, Vr) { return Math.abs(d) < 1e-12 ? 0 : Math.sqrt(4 * Math.abs(d) * Pu / Vr); };
    var Fr = Math.sqrt(Fo * Fo + Fi * Fi);
    return { Pu: Pu, Rlat: Rl, Rvert: Rv, Fout: Fo, Fin: Fi, Fres: Fr, dcLat: dcl, dcVert: dcv,
             VrLat: Vl, VrVert: Vv, okLat: Fo <= Vl + 1e-9, okVert: Fi <= Vv + 1e-9,
             okRes: Fr <= Math.min(Vl, Vv) + 1e-9, RminOk: Math.min(Rl, Rv) >= 6000 - 1e-9,
             LtMinLat: ltm(dx, Vl), LtMinVert: ltm(dy, Vv), sMax: Math.min(3 * od, 600) };
  };

  // 相鄰（上下疊放）曲線管道互推（同 tendon_profile.adjacent_duct_radial_check）：AASHTO 5.10.4.3.1
  // 三擇一——加大間距使 V_n 足夠／配圍束筋承擔徑向力／內側管先灌漿再拉外側。條文為**定性**，
  // 故：路徑一以 Eq.3 估（d_c 取「淨間距+½管徑」，條文未定義管對管 d_c，屬解釋）；
  // 路徑二依 §5.10.4.3 服務狀態 f_s ≤ 0.6f_y（f_y ≤ 420）回推鋼筋量。疊放時上方承受該列之和。
  BC.adjacentDuctRadialCheck = function (Pu, Ps, R, nStack, sClear, fci, od, fy, phi, forceToward) {
    od = od == null ? 100 : od; fy = Math.min(fy == null ? 420 : fy, 420);
    phi = phi == null ? 0.9 : phi; forceToward = forceToward !== false;
    var Fe = (R && R !== Infinity) ? Pu / R : 0, dc = sClear + od / 2,
        Vr = phi * 0.33 * Math.sqrt(fci) * dc,
        sReq = Math.max(0, Fe / (phi * 0.33 * Math.sqrt(fci)) - od / 2),
        FsStack = ((R && R !== Infinity) ? Ps / R : 0) * nStack,
        AsPerMm = FsStack / (0.6 * fy), sMax = Math.min(3 * od, 600);
    return { applies: !!(forceToward && nStack >= 2 && Fe > 0), nStack: nStack, Feach: Fe,
             Fstack: Fe * nStack, sClear: sClear, dcBetween: dc, VrBetween: Vr,
             spacingOk: Fe <= Vr + 1e-9, sReq: sReq, AsPerMm: AsPerMm, tieSmax: sMax,
             AsPerTie: AsPerMm * sMax, fyUsed: fy };
  };
  // 台灣 §8.9.1／§8.9.2 箱梁翼板最小厚度（同 model.box_slab_thickness_TW）
  // ⚠ RC 箱梁底板才是「梁腹間淨距/16」（第七章 §7.1.22 10.(3)b）；PC 頂底板同為 1/30，勿混用。
  BC.boxSlabThicknessTW = function (topT, botT, clearSpan, precast) {
    var tf = precast ? 140 : 150, bf = precast ? 130 : 140,
        tr = Math.max(clearSpan / 30, tf), br = Math.max(clearSpan / 30, bf);
    return { clearSpan: clearSpan, topReq: tr, botReq: br, topOk: topT >= tr - 1e-9,
             botOk: botT >= br - 1e-9, rcBotReq: Math.max(clearSpan / 16, 140) };
  };

  // ── 腹板抗剪 D1 ───────────────────────────────────────
  BC.principalTensionLimitTW = function (fc) { return 0.094 * sqrt(fc); };
  BC.shearWeb = function (Pe, sec, e, fc, bw, dv, Vu, xCtrl, L, nWebs, phi, fsy) {
    return BC.shearWebAt(Pe, sec, (4 * e / L) * (1 - 2 * xCtrl / L), fc, bw, dv, Vu, nWebs, phi, fsy);
  };
  // 任意斷面（連續梁墩兩側等）：slope＝de/dx（帶號）、Vu 帶號；Vp 僅當 slope 與 Vu 同號時有利，否則為負（同 shear.shear_web_at）
  BC.shearWebAt = function (Pe, sec, slope, fc, bw, dv, Vu, nWebs, phi, fsy) {
    nWebs = nWebs || 2; phi = phi || 0.85; fsy = fsy || 420;
    var fpc = Pe / sec.A, Va = Math.abs(Vu), Vp = (Pe / nWebs) * Math.abs(slope) * (slope * Vu >= 0 ? 1 : -1);
    var tau = Va / (bw * dv), sigma1 = -fpc / 2 + sqrt(Math.pow(fpc / 2, 2) + tau * tau);
    var lim = BC.principalTensionLimitTW(fc), vc = 0.094 * sqrt(fc);
    var Vcw = (vc + 0.3 * fpc) * bw * dv + Vp, Vs_req = (Va - phi * Vcw) / phi;
    return { fpc: fpc, slope: slope, Vp: Vp, tau: tau, sigma1: sigma1, sigma1_limit: lim,
             sigma1_ok: sigma1 <= lim, Vcw: Vcw, Vs_req: Vs_req, Av_s_req: max(Vs_req, 0) / (fsy * dv) };
  };
  BC.phiVn = function (Vcw, Av_s, dv, phi, fsy) { phi = phi || 0.85; fsy = fsy || 420; return phi * (Vcw + Av_s * fsy * dv); };
  // 台灣最小腹板鋼筋 Av/s = 0.345·b′/fsy（§8.20.3 3. 式 8-31，原文 kgf 制 3.5；第七章 §7.1.9 1.(2) 同式）。
  // 🔴 2026-09-20 更正：原寫 max(0.2√f'c·b/fsy, 0.35·b/fsy) 有兩個錯——√f'c 項非台灣規範（ACI／AASHTO），
  //    且 0.2 是 kgf/cm² 制係數卻餵 MPa（ACI SI 0.062／AASHTO 0.083），f'c=40 時為明文值的 3.7 倍。
  //    使用者裁示：判定用台灣明文，AASHTO 5.8.2.5 以 AvSminAASHTO 並列參考。fc 參數保留相容。
  BC.AvSminTW = function (fc, bw, fsy) { fsy = fsy || 420; return 0.345 * bw / fsy; };
  BC.AvSminAASHTO = function (fc, bw, fy) { fy = fy || 420; return 0.083 * sqrt(fc) * bw / fy; };

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

  // ── 台灣第十二章耐久性保護層（同 durability.durability_cover；原文 p.341–343）────
  // 表 12.1：等級 I＝箱梁內部；II＝乾濕交替（外露面）。取 max(表列, §8.25.1 基本 40 mm)。
  BC.TW_COVER_GENERAL = { 50: { I: [[0.50, 30], [0.45, 25]], II: [[0.50, 40], [0.45, 35], [0.40, 30]] },
                          100: { I: [[0.50, 35], [0.45, 30]], II: [[0.45, 45], [0.40, 40], [0.35, 35]] } };
  BC.TW_SALT_MAX_WC = { '極嚴重': 0.40, '嚴重': 0.40, '中度': 0.45 };
  BC.TW_COVER_SALT = {
    '基礎基樁': { 50: [100, 100, 100], 100: [100, 100, 100] }, '柱牆': { 50: [100, 75, 75], 100: [100, 100, 75] },
    '橋面版頂層筋': { 50: [65, 55, 50], 100: [75, 65, 60] }, '橋面版下層筋': { 50: [65, 55, 50], 100: [75, 65, 60] },
    '箱梁底層筋': { 50: [65, 55, 50], 100: [75, 65, 60] }, '梁腹版外露面': { 50: [65, 55, 50], 100: [75, 65, 60] },
    '未曝露面': { 50: [40, 40, 40], 100: [40, 40, 40] } };
  BC.durabilityCover = function (env, life, grade, wc, member, zone, basic) {
    env = env || 'general'; life = life || 50; grade = grade || 'II'; wc = wc == null ? 0.45 : wc;
    member = member || '梁腹版外露面'; zone = zone || '中度'; basic = basic == null ? 40 : basic;
    var tbl = null, src;
    if (env === 'general') {
      var ok = BC.TW_COVER_GENERAL[life][grade].filter(function (r) { return wc <= r[0] + 1e-9; }).map(function (r) { return r[1]; });
      tbl = ok.length ? Math.min.apply(null, ok) : null; src = '表12.2 一般環境 等級' + grade + ' ' + life + '年';
    } else {
      if (wc <= BC.TW_SALT_MAX_WC[zone] + 1e-9) tbl = BC.TW_COVER_SALT[member][life][{ '極嚴重': 0, '嚴重': 1, '中度': 2 }[zone]];
      src = '表12.5 鹽害 ' + member + ' ' + zone + '鹽害區 ' + life + '年';
    }
    return { table: tbl, basic: basic, required: tbl == null ? null : max(tbl, basic), wc_ok: tbl != null, source: src };
  };

  // ── 套管尺寸檢核（同 tendon_profile.duct_size_check）─────────────
  // 表 8.3 與 PTI Table 4.4 皆為**內徑**；外徑依廠商（PTI §4.4.5）。排列用外徑、面積比用內徑。
  BC.TW_DUCT_MAX_ID = { '12.7': { 22: 90, 19: 90, 12: 75, 7: 55 }, '15.2': { 22: 110, 19: 100, 12: 85, 7: 70 } };
  BC.DUCT_AREA_RATIO = { tw: 2.0, pti_push: 2.25, pti_pull: 2.5, pti_short: 2.0 };
  BC.ductSizeCheck = function (nStrands, ductId, ductOd, strandDia, strandArea, rule) {
    strandDia = strandDia || 15.2; strandArea = strandArea || 140; rule = rule || 'tw';
    var Aps = nStrands * strandArea, Ad = Math.PI * ductId * ductId / 4, req = BC.DUCT_AREA_RATIO[rule],
        tb = BC.TW_DUCT_MAX_ID[strandDia.toFixed(1)], idm = tb && tb[nStrands] != null ? tb[nStrands] : null,
        od = ductOd == null ? ductId : ductOd;
    return { A_ps: Aps, A_duct: Ad, ratio: Ad / Aps, ratio_req: req, area_ok: Ad / Aps >= req - 1e-9,
             id_max_tw: idm, id_ok: idm == null || ductId <= idm + 1e-9,
             od_gt_id: od > ductId + 1e-9, wall: (od - ductId) / 2 };
  };

  // ── 中間錨碇齒塊（blister）錨碇區（同 bridgecalc.blister）───────
  // 齒塊 ≠ 端部錨碇 ≠ 轉向塊：後方無鋼腱延伸提供回力 → Tie-back 是唯一抵抗
  // 後向分離的機制；介面剪力傳的是 P·cosα（近全腱力）而非轉向塊的 P·sinα。
  BC.blisterLocalBearing = function (Pd, A_plate, A_bearing, fci, phiB) {
    phiB = phiB == null ? 0.70 : phiB;
    var fb = Pd * 1e3 / A_plate, raw = 0.85 * phiB * fci * sqrt(A_bearing / A_plate),
        cap = 1.5 * fci, allow = min(raw, cap);
    return { f_b: fb, f_b_allow: allow, f_b_uncapped: raw, capped: raw > cap,
             ratio: fb / allow, ok: fb <= allow, spiral_factor_req: max(1, fb / allow) };
  };
  BC.blisterBursting = function (Pu, a_plate, h, e, fy, phi) {
    e = e || 0; fy = fy || 420; phi = phi == null ? 0.9 : phi;
    var F = 0.25 * Pu * (1 - a_plate / h), d = 0.5 * (h - 2 * abs(e));
    return { F_burst: F, d_burst: d, As_burst: F * 1e3 / (phi * fy),
             zone_from: 0.5 * d, zone_to: 1.5 * d };
  };
  BC.blisterTieback = function (Ps, f_cb, A_cb, fy, a_plate) {
    f_cb = f_cb || 0; A_cb = A_cb || 0; fy = fy || 420; a_plate = a_plate || 200;
    var T = 0.25 * Ps, C = f_cb * A_cb / 1e3, fs = min(0.6 * fy, 248);
    return { T_req: T, C_precomp: C, fs_allow: fs, As: max(0, T - C) * 1e3 / fs,
             As_conservative: T * 1e3 / fs, max_dist_from_axis: a_plate };
  };
  BC.blisterSpalling = function (Ps, fy, phi, ratio) {
    fy = fy || 420; phi = phi == null ? 0.9 : phi; ratio = ratio == null ? 0.02 : ratio;
    var F = ratio * Ps;
    return { F_spall: F, As_spall: F * 1e3 / (phi * fy) };
  };
  // ⚠ 剪力摩擦有**面積上限**：As_vf 算得出來不代表做得到。介面面積不足時
  //   加多少鋼筋都沒用，只能加大齒塊。K1/K2 隨規範版次不同。
  BC.blisterInterfaceShear = function (Ps, alphaDeg, mu, fsd, Nd, gamma0, A_int, fc, K1, K2) {
    alphaDeg = alphaDeg == null ? 5 : alphaDeg; mu = mu || 1.0; fsd = fsd || 360;
    Nd = Nd || 0; gamma0 = gamma0 || 1.0; A_int = A_int || 0; fc = fc || 40;
    K1 = K1 == null ? 0.25 : K1; K2 = K2 == null ? 10.3 : K2;
    var V = Ps * Math.cos(alphaDeg * Math.PI / 180),
        As = max(0, gamma0 * V * 1e3 / mu + Nd * 1e3) / fsd,
        tauCap = min(K1 * fc, K2), Vcap = A_int ? tauCap * A_int / 1e3 : 0;
    return { V_int: V, As_vf: As, tau: A_int ? V * 1e3 / A_int : 0,
             tau_cap: tauCap, V_cap: Vcap, area_ok: A_int > 0 && V <= Vcap + 1e-9,
             A_req: tauCap > 0 ? V * 1e3 / tauCap : Infinity };
  };
  BC.blisterGeometryCheck = function (L_b, W_b, a_plate, b_plate, straightHave, straightReq, edgeReq) {
    straightReq = straightReq == null ? 400 : straightReq;
    edgeReq = edgeReq == null ? 50 : edgeReq;
    var Lmin = a_plate + 2 * edgeReq, edgeHave = (W_b - b_plate) / 2,
        sOk = straightHave >= straightReq - 1e-9, eOk = edgeHave >= edgeReq - 1e-9,
        lOk = L_b >= Lmin - 1e-9;
    return { straight_req: straightReq, straight_have: straightHave, straight_ok: sOk,
             edge_req: edgeReq, edge_have: edgeHave, edge_ok: eOk,
             L_min: Lmin, L_ok: lOk, ok: sOk && eOk && lOk };
  };
  // 三個設計力刻意分開：P_s 服務(Tie-back/剝裂/介面)、P_u 係數化(爆裂)、P_d 張拉(局部承壓)
  BC.blisterDesign = function (Ps, Pu, Pd, o) {
    o = o || {};
    var ap = o.a_plate == null ? 200 : o.a_plate, bp = o.b_plate == null ? 200 : o.b_plate,
        Lb = o.L_b == null ? 400 : o.L_b, Wb = o.W_b == null ? 300 : o.W_b,
        Db = o.D_b == null ? 250 : o.D_b, Ab = o.A_bearing == null ? 60000 : o.A_bearing,
        fci = o.fci == null ? 35 : o.fci, fcb = o.f_cb || 0, Acb = o.A_cb || 0,
        al = o.alpha_deg == null ? 5 : o.alpha_deg, mu = o.mu || 1.0,
        fy = o.fy || 420, fsd = o.fsd || 360, fc = o.fc || 40,
        sh = o.straight_have == null ? 400 : o.straight_have;
    var Aplate = ap * bp, Aint = Lb * Db;
    var bearing = BC.blisterLocalBearing(Pd, Aplate, Ab, fci),
        burst = BC.blisterBursting(Pu, ap, Lb, 0, fy),
        tie = BC.blisterTieback(Ps, fcb, Acb, fy, ap),
        spall = BC.blisterSpalling(Ps, fy),
        face = BC.blisterInterfaceShear(Ps, al, mu, fsd, 0, 1, Aint, fc),
        geom = BC.blisterGeometryCheck(Lb, Wb, ap, bp, sh);
    var items = { '爆裂': burst.As_burst, 'Tie-back': tie.As_conservative,
                  '剝裂': spall.As_spall, '介面剪力摩擦': face.As_vf };
    var gov = Object.keys(items).reduce(function (a, b) { return items[b] > items[a] ? b : a; });
    var tot = Object.keys(items).reduce(function (a, k) { return a + items[k]; }, 0);
    return { bearing: bearing, burst: burst, tie: tie, spall: spall, face: face,
             geom: geom, As_total_conservative: tot, governing: gov, A_interface: Aint };
  };

  // ── STM 壓桿／節點（同 bridgecalc.stm）：AASHTO LRFD 2008 SI 原文式 ──────────
  // 🔴 AASHTO 2008 **沒有 β_s 表**：f_cu = f'c/(0.8+170ε₁) ≤ 0.85f'c（Eq. 5.6.3.3.3-1）、
  //    ε₁ = ε_s+(ε_s+0.002)cot²α_s（-2）；節點 CCC 0.85／CCT 0.75／CTT 0.65 ×φf'c（5.6.3.5）、φ=0.70。
  //    ACI 318 式 β 表（節點 1.00/0.80/0.60/0.40、壓桿 0.60/0.75/1.00）另列對照，兩者僅 CCC 同值。
  BC.NODE_LIMIT_AASHTO = { CCC: 0.85, CCT: 0.75, CTT: 0.65 };
  BC.BETA_NODE_ACI = { CCC: 1.00, CCT: 0.80, CTT: 0.60, TTT: 0.40 };
  BC.BETA_STRUT_ACI = { unconfined: 0.60, reinforced: 0.75, confined: 1.00 };
  BC.EPS_S_YIELD = 420 / 200000;
  BC.strutEps1 = function (epsS, alphaDeg) {
    var a = alphaDeg * Math.PI / 180;
    return epsS + (epsS + 0.002) / (Math.tan(a) * Math.tan(a));
  };
  BC.strutFcuAASHTO = function (fc, epsS, alphaDeg) {
    epsS = epsS == null ? BC.EPS_S_YIELD : epsS; alphaDeg = alphaDeg == null ? 45 : alphaDeg;
    return Math.min(fc / (0.8 + 170 * BC.strutEps1(epsS, alphaDeg)), 0.85 * fc);
  };
  BC.nodeCapacityAASHTO = function (fc, type, An, phi) {
    return BC.NODE_LIMIT_AASHTO[type] * (phi == null ? 0.70 : phi) * fc * An;
  };
  BC.fCuACI = function (fc, beta) { return 0.85 * beta * fc; };
  BC.nodeCapacityACI = function (fc, type, An, phi) {
    return (phi == null ? 0.70 : phi) * BC.fCuACI(fc, BC.BETA_NODE_ACI[type]) * An;
  };
  // 齒塊 STM（同 blister.blister_stm）：壓桿寬 w_s = l_b·sinθ + w_t·cosθ（CCT 節點幾何）。
  // ⚠ 錨板面屬局部區（5.10.9.7，見 blisterLocalBearing）；本檢核針對齒塊**根部節點**，A_node 請給實際面積。
  BC.blisterSTM = function (C_kN, Pnode_kN, o) {
    o = o || {};
    var ap = o.a_plate == null ? 200 : o.a_plate, bp = o.b_plate == null ? 200 : o.b_plate,
        wt = o.w_tie == null ? 150 : o.w_tie, th = o.theta_deg == null ? 45 : o.theta_deg,
        webT = o.web_t == null ? 350 : o.web_t, fc = o.fc == null ? 40 : o.fc,
        nt = o.node_type || 'CCT', es = o.eps_s == null ? BC.EPS_S_YIELD : o.eps_s,
        bs = o.beta_s_aci || 'reinforced', phi = o.phi == null ? 0.70 : o.phi;
    var r = th * Math.PI / 180, ws = ap * Math.sin(r) + wt * Math.cos(r),
        bEff = Math.min(bp, webT), Acs = ws * bEff, C = C_kN * 1e3, Pn = Pnode_kN * 1e3;
    var fcuA = BC.strutFcuAASHTO(fc, es, th), fcuI = BC.fCuACI(fc, BC.BETA_STRUT_ACI[bs]);
    var An = o.A_node == null ? ap * bEff : o.A_node;
    var FnnA = BC.nodeCapacityAASHTO(fc, nt, An, phi), FnnI = BC.nodeCapacityACI(fc, nt, An, phi);
    return { C_strut: C, theta_deg: th, w_s: ws, b_eff: bEff, A_cs: Acs,
             fcu_aashto: fcuA, phiFns_aashto: phi * fcuA * Acs, strut_ok_aashto: phi * fcuA * Acs >= C,
             fcu_aci: fcuI, phiFns_aci: phi * fcuI * Acs, strut_ok_aci: phi * fcuI * Acs >= C,
             A_n: An, A_n_req: Pn / (BC.NODE_LIMIT_AASHTO[nt] * phi * fc), A_cs_req: C / (phi * fcuA),
             node_type: nt, phiFnn_aashto: FnnA, node_ok_aashto: FnnA >= Pn,
             phiFnn_aci: FnnI, node_ok_aci: FnnI >= Pn, eps1: BC.strutEps1(es, th),
             util_strut: C / (phi * fcuA * Acs), util_node: Pn / FnnA };
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

  // ── 連續梁變斷面（中墩底板加厚）：變 EI 柔度（同 variable_section.py）──
  // F_ij＝∫m_i·m_j/I_rel；點載重轉角項 δ_i(p)＝u_i(p)：簡支跨以曲率 κ_i＝m_i/I_rel 的撓度（Maxwell 互易）
  var GL3 = [[-Math.sqrt(3 / 5), 5 / 9], [0, 8 / 9], [Math.sqrt(3 / 5), 5 / 9]];
  function contCells(a, b, cuts, h) {
    var seen = {}, pts = [];
    [a, b].concat(cuts.filter(function (c) { return c > a && c < b; })).forEach(function (v) { if (!seen[v]) { seen[v] = 1; pts.push(v); } });
    pts.sort(function (x, y) { return x - y; });
    var out = [];
    for (var i = 0; i < pts.length - 1; i++) {
      var u = pts[i], v = pts[i + 1], n = Math.max(1, Math.round((v - u) / h)), d = (v - u) / n;
      for (var k = 0; k < n; k++) out.push([u + k * d, u + (k + 1) * d]);
    }
    return out;
  }
  BC.gaussIntegrate = function (f, a, b, cuts, h) {
    var s = 0;
    contCells(a, b, cuts || [], h || 0.25).forEach(function (c) {
      var m = 0.5 * (c[0] + c[1]), r = 0.5 * (c[1] - c[0]);
      s += r * (GL3[0][1] * f(m + r * GL3[0][0]) + GL3[1][1] * f(m + r * GL3[1][0]) + GL3[2][1] * f(m + r * GL3[2][0]));
    });
    return s;
  };
  BC.haunchProfile = function (spans, topW, topT, botW, botT, webT, nWeb, h, botTPier, length) {
    var xs = [0], lens = [], brk = [], seen = {}, j;
    spans.forEach(function (L) { xs.push(xs[xs.length - 1] + L); });
    for (j = 1; j < spans.length; j++) {
      var b = Math.max(0, Math.min(length, 0.5 * spans[j - 1], 0.5 * spans[j]));
      lens.push(b); [xs[j] - b, xs[j], xs[j] + b].forEach(function (v) { v = Math.round(v * 1e9) / 1e9; if (!seen[v]) { seen[v] = 1; brk.push(v); } });
    }
    brk.sort(function (a, c) { return a - c; });
    var ref = BC.sectionFromDims(topW, topT, botW, botT, webT, nWeb, h), tPier = Math.max(botT, botTPier), cache = {};
    var P = { spans: spans.slice(), xs: xs, lengths: lens, breaks: brk, secRef: ref, botTPier: tPier };
    P.botTAt = function (x) {
      var t = botT;
      for (var i = 1; i < xs.length - 1; i++) { var d = Math.abs(x - xs[i]), bb = lens[i - 1];
        if (bb > 0 && d < bb) t = Math.max(t, tPier + (botT - tPier) * d / bb); }
      return t;
    };
    P.sectionAt = function (x) { var k = Math.round(x * 1e6); return cache[k] || (cache[k] = BC.sectionFromDims(topW, topT, botW, P.botTAt(x), webT, nWeb, h)); };
    P.Irel = function (x) { return P.sectionAt(x).I / ref.I; };
    P.Arel = function (x) { return P.sectionAt(x).A / ref.A; };
    P.dyb = function (x) { return ref.yb - P.sectionAt(x).yb; };
    return P;
  };
  BC.contFlex = function (spans, Irel, breaks, grid) {
    breaks = breaks || []; grid = grid || 0.05;
    var xs = [0], n = spans.length, nu = n - 1, i, j;
    spans.forEach(function (L) { xs.push(xs[xs.length - 1] + L); });
    function mh(i, x) { var a = xs[i - 1], c = xs[i], b = xs[i + 1];
      if (x >= a - 1e-12 && x <= c) return (x - a) / (c - a); if (x >= c && x <= b + 1e-12) return (b - x) / (b - c); return 0; }
    var F = [], cuts = xs.concat(breaks);
    for (i = 0; i < nu; i++) { F.push([]); for (j = 0; j < nu; j++) F[i].push(0); }
    for (i = 0; i < nu; i++) {
      (function (i) {
        F[i][i] = BC.gaussIntegrate(function (x) { var m = mh(i + 1, x); return m * m / Irel(x); }, xs[i], xs[i + 2], cuts, 0.25);
        if (i + 1 < nu) F[i][i + 1] = F[i + 1][i] = BC.gaussIntegrate(function (x) { return mh(i + 1, x) * mh(i + 2, x) / Irel(x); }, xs[i + 1], xs[i + 2], cuts, 0.25);
      })(i);
    }
    var tab = {};
    function table(i, k) {
      var x0 = xs[k], L = spans[k], cells = contCells(x0, x0 + L, breaks, grid), ts = [x0], S1 = [0], S2 = [0];
      cells.forEach(function (c) {
        var u = c[0], v = c[1], m = 0.5 * (u + v), r = 0.5 * (v - u), k1 = 0, k2 = 0;
        GL3.forEach(function (gw) { var t = m + r * gw[0], kap = mh(i, t) / Irel(t); k1 += r * gw[1] * kap; k2 += r * gw[1] * kap * (v - t); });
        S2.push(S2[S2.length - 1] + (v - u) * S1[S1.length - 1] + k2); S1.push(S1[S1.length - 1] + k1); ts.push(v);
      });
      var cc = S2[S2.length - 1] / L;
      return { ts: ts, u: ts.map(function (t, q) { return -S2[q] + (t - x0) * cc; }), du: S1.map(function (s1) { return -s1 + cc; }) };
    }
    for (i = 1; i < n; i++) { tab[i + ',' + (i - 1)] = table(i, i - 1); tab[i + ',' + i] = table(i, i); }
    function uAt(i, k, p) {
      var T = tab[i + ',' + k], ts = T.ts, lo = 0, hi = ts.length - 1;
      if (p <= ts[0]) return T.u[0]; if (p >= ts[hi]) return T.u[hi];
      while (hi - lo > 1) { var md = (lo + hi) >> 1; if (ts[md] <= p) lo = md; else hi = md; }
      var hh = ts[hi] - ts[lo], s = (p - ts[lo]) / hh, s2 = s * s, s3 = s2 * s;
      return (2 * s3 - 3 * s2 + 1) * T.u[lo] + (s3 - 2 * s2 + s) * hh * T.du[lo] + (-2 * s3 + 3 * s2) * T.u[hi] + (s3 - s2) * hh * T.du[hi];
    }
    function solve(delta) {
      if (nu <= 0) { var z = []; for (var q = 0; q <= n; q++) z.push(0); return z; }
      var A = F.map(function (row, r) { return row.slice().concat([-delta[r]]); }), X = [], a, b2, c;
      for (a = 0; a < nu; a++) for (b2 = a + 1; b2 < nu; b2++) { var f = A[b2][a] / A[a][a]; for (c = a; c <= nu; c++) A[b2][c] -= f * A[a][c]; }
      for (a = 0; a < nu; a++) X.push(0);
      for (a = nu - 1; a >= 0; a--) { var sm = A[a][nu]; for (c = a + 1; c < nu; c++) sm -= A[a][c] * X[c]; X[a] = sm / A[a][a]; }
      return [0].concat(X, [0]);
    }
    function spanOf(x) { for (var q = 0; q < n; q++) if (x <= xs[q + 1] + 1e-9) return q; return n - 1; }
    var o = { F: F, xs: xs, spans: spans, breaks: breaks };
    o.supportMomentsPoint = function (p, P) {
      P = P == null ? 1 : P;
      if (n < 2 || p < -1e-9 || p > xs[n] + 1e-9) { var z = []; for (var q = 0; q <= n; q++) z.push(0); return z; }
      var k = spanOf(p), delta = []; for (var q2 = 0; q2 < nu; q2++) delta.push(0);
      [k, k + 1].forEach(function (ii) { if (ii >= 1 && ii <= n - 1) delta[ii - 1] = P * uAt(ii, k, p); });
      return solve(delta);
    };
    o.supportMomentsDist = function (w) {
      var delta = []; for (var q = 0; q < nu; q++) delta.push(0);
      for (var ii = 1; ii < n; ii++) [ii - 1, ii].forEach(function (k) {
        delta[ii - 1] += BC.gaussIntegrate(function (t) { return w(t) * uAt(ii, k, t); }, xs[k], xs[k + 1], breaks, 0.25); });
      return solve(delta);
    };
    o.momentAt = function (X, x, M0) { var q = spanOf(x), t = (x - xs[q]) / spans[q]; return M0 + X[q] * (1 - t) + X[q + 1] * t; };
    o.M0Dist = function (w, x) {
      var k = spanOf(x), x0 = xs[k], L = spans[k], a = x - x0;
      return BC.gaussIntegrate(function (t) { var s2 = t - x0; return (s2 <= a ? s2 * (L - a) / L : a * (L - s2) / L) * w(t); }, x0, x0 + L, breaks.concat([x]), 0.25);
    };
    o.spanOf = spanOf;
    return o;
  };

  // ── 連續梁影響線與台灣 HS20-44 彎矩包絡（解析三彎矩法，EI 常數；同 influence_cont.py）──
  // 點載重 P 距跨左端 a（b=L−a）：∫M0·m＝P·a·b·(L+a)/(6L)（右端支承）、P·a·b·(L+b)/(6L)（左端）
  // 台灣 §3.9：負彎矩車道另加一個等集中載重於他跨（共 2 個）；§3.13 衝擊 L：正彎矩＝該跨、負彎矩＝相鄰兩跨平均
  function contXs(spans) { var xs = [0]; spans.forEach(function (L) { xs.push(xs[xs.length - 1] + L); }); return xs; }
  function contSpanOf(xs, x) { var n = xs.length - 1; for (var i = 0; i < n; i++) if (x <= xs[i + 1] + 1e-9) return i; return n - 1; }
  function contTridiag(spans, rhs) {
    var n = spans.length, nu = n - 1, i;
    if (nu <= 0) { var z = []; for (i = 0; i <= n; i++) z.push(0); return z; }
    var a = [], b = [], c = [], d = rhs.slice();
    for (i = 0; i < nu; i++) { a.push(spans[i]); b.push(2 * (spans[i] + spans[i + 1])); c.push(spans[i + 1]); }
    for (i = 1; i < nu; i++) { var w = a[i] / b[i - 1]; b[i] -= w * c[i - 1]; d[i] -= w * d[i - 1]; }
    var X = new Array(nu); X[nu - 1] = d[nu - 1] / b[nu - 1];
    for (i = nu - 2; i >= 0; i--) X[i] = (d[i] - c[i] * X[i + 1]) / b[i];
    return [0].concat(X, [0]);
  }
  BC.contSupportMomentsPoint = function (spans, p, P) {
    P = P == null ? 1 : P;
    var xs = contXs(spans), n = spans.length, rhs = [], i;
    for (i = 0; i < n - 1; i++) rhs.push(0);
    if (n < 2 || p < -1e-9 || p > xs[n] + 1e-9) { var z = []; for (i = 0; i <= n; i++) z.push(0); return z; }
    var k = contSpanOf(xs, p), L = spans[k], aa = Math.min(Math.max(p - xs[k], 0), L), bb = L - aa;
    if (k >= 1) rhs[k - 1] += -6 * P * aa * bb * (L + bb) / (6 * L);
    if (k <= n - 2) rhs[k] += -6 * P * aa * bb * (L + aa) / (6 * L);
    return contTridiag(spans, rhs);
  };
  BC.contSupportMomentsUniform = function (spans, w) {
    var n = spans.length, rhs = [], k;
    for (k = 0; k < n - 1; k++) rhs.push(0);
    for (k = 0; k < n; k++) { var t = -6 * w[k] * Math.pow(spans[k], 3) / 24;
      if (k >= 1) rhs[k - 1] += t; if (k <= n - 2) rhs[k] += t; }
    return contTridiag(spans, rhs);
  };
  function contMomentAt(spans, xs, X, x, M0) { var i = contSpanOf(xs, x), t = (x - xs[i]) / spans[i]; return M0 + X[i] * (1 - t) + X[i + 1] * t; }
  function contEta(spans, xs, X, x, p) {
    var i = contSpanOf(xs, x), k = contSpanOf(xs, p), M0 = 0;
    if (i === k) { var L = spans[i], a = x - xs[i], q = p - xs[i]; M0 = q <= a ? q * (L - a) / L : a * (L - q) / L; }
    return contMomentAt(spans, xs, X, x, M0);
  }
  // 強度組合恆載因數依效應正負號取 max／min（AASHTO γ_p；同 influence_cont.factored）。
  // 恆載與所求效應同號＝不利 → γ_max；反號＝有利 → γ_min。控制斷面不受影響。
  BC.G_DC = 1.25; BC.G_DC_MIN = 0.90; BC.G_DW = 1.50; BC.G_DW_MIN = 0.65; BC.G_LL = 1.75;
  BC.factored = function (dc, dw, ll, wantPos, gdc, gdcm, gdw, gdwm, gll) {
    gdc = gdc == null ? BC.G_DC : gdc; gdcm = gdcm == null ? BC.G_DC_MIN : gdcm;
    gdw = gdw == null ? BC.G_DW : gdw; gdwm = gdwm == null ? BC.G_DW_MIN : gdwm; gll = gll == null ? BC.G_LL : gll;
    var f = function (v, gmax, gmin) { return ((v >= 0) === wantPos ? gmax : gmin) * v; };
    return f(dc, gdc, gdcm) + f(dw, gdw, gdwm) + gll * ll;
  };
  // 部分均布載重（跨內任意區段）：三彎矩右端項由集中載重式對載重長度解析積分
  //   左端 −(w/L)[L²a² − L a³ + a⁴/4]、右端 −(w/L)[L²a²/2 − a⁴/4]（a1=0、a2=L 時各得 −wL³/4）
  BC.contSupportMomentsPartial = function (spans, loads) {
    var xs = contXs(spans), n = spans.length, rhs = [], i, z;
    if (n < 2) { z = []; for (i = 0; i <= n; i++) z.push(0); return z; }
    for (i = 0; i < n - 1; i++) rhs.push(0);
    loads.forEach(function (ld) {
      var x1 = ld[0], x2 = ld[1], w = ld[2];
      if (w === 0 || x2 <= x1) return;
      for (var k = 0; k < n; k++) {
        var L = spans[k], a1 = min(max(x1 - xs[k], 0), L), a2 = min(max(x2 - xs[k], 0), L);
        if (a2 <= a1) continue;
        var tl = -(w / L) * ((L * L * a2 * a2 - L * a2 * a2 * a2 + Math.pow(a2, 4) / 4)
                             - (L * L * a1 * a1 - L * a1 * a1 * a1 + Math.pow(a1, 4) / 4));
        var tr = -(w / L) * ((L * L * a2 * a2 / 2 - Math.pow(a2, 4) / 4)
                             - (L * L * a1 * a1 / 2 - Math.pow(a1, 4) / 4));
        if (k >= 1) rhs[k - 1] += tl;
        if (k <= n - 2) rhs[k] += tr;
      }
    });
    return contTridiag(spans, rhs);
  };
  function M0Partial(L, a, local) {
    var m = 0;
    local.forEach(function (ld) {
      var a1 = ld[0], a2 = ld[1], w = ld[2];
      if (a2 <= a1 || w === 0) return;
      var W = w * (a2 - a1), cg = (a1 + a2) / 2, R = W * (L - cg) / L;
      if (a <= a1) m += R * a;
      else if (a >= a2) m += R * a - W * (a - cg);
      else m += R * a - w * (a - a1) * (a - a1) / 2;
    });
    return m;
  }
  BC.contMomentPartialUDL = function (spans, loads, x) {
    var xs = contXs(spans), X = BC.contSupportMomentsPartial(spans, loads),
        i = contSpanOf(xs, x), L = spans[i], a = x - xs[i], local = [];
    loads.forEach(function (ld) {
      var a1 = min(max(ld[0] - xs[i], 0), L), a2 = min(max(ld[1] - xs[i], 0), L);
      if (a2 > a1) local.push([a1, a2, ld[2]]);
    });
    return contMomentAt(spans, xs, X, x, M0Partial(L, a, local));
  };
  BC.contMomentIL = function (spans, x, p) { var xs = contXs(spans); return contEta(spans, xs, BC.contSupportMomentsPoint(spans, p), x, p); };
  BC.contDLMoment = function (spans, w, x) {
    var xs = contXs(spans), X = BC.contSupportMomentsUniform(spans, spans.map(function () { return w; }));
    var i = contSpanOf(xs, x), a = x - xs[i];
    return contMomentAt(spans, xs, X, x, w * a * (spans[i] - a) / 2);
  };
  // ── 施工階段體系轉換的潛變重分配（同 bridgecalc.staging）──────
  // 載重在施工時作用於體系 I（逐跨施工＝各跨簡支），連續後潛變使內力朝
  // 「一開始就作用在體系 II」的解漂移：M(∞) = M_I + λ(M_II − M_I)。
  // ⚠ 決定 λ 的是**剩餘**潛變 Δφ = φ(∞) − φ(t₁)，故越早合龍重分配越多。
  // phiR（嚴格 AAEM）：分母用束制彎矩自身的潛變 φ(∞,t1)；不給時以 Δφ 近似（λ 偏大、偏保守）
  BC.redistributionFactor = function (dphi, chi, method, phiR) {
    chi = chi == null ? 0.8 : chi; method = method || 'trost';
    if (dphi < 0) throw new Error('Δφ 不可為負');
    var pr = phiR == null ? dphi : phiR;
    var lam = method === 'dischinger' ? 1 - exp(-dphi) : dphi / (1 + chi * pr);
    var valid = lam >= 0 && lam <= 1 + 1e-12;
    return { lam: lam, dphi: dphi, chi: chi, method: method, valid: valid,
             note: valid ? '' : 'λ=' + lam.toFixed(3) + ' 超出 [0,1]，Δφ 過大已逸出模型適用範圍' };
  };
  BC.creepRedistribution = function (xs, M_I, M_II, dphi, chi, method, phiR) {
    if (!(xs.length === M_I.length && xs.length === M_II.length)) throw new Error('長度須相同');
    var f = BC.redistributionFactor(dphi, chi, method, phiR), pts = [];
    for (var i = 0; i < xs.length; i++) {
      var d = f.lam * (M_II[i] - M_I[i]);
      pts.push({ x: xs[i], M_I: M_I[i], M_II: M_II[i], dM: d, M_inf: M_I[i] + d });
    }
    var hi = pts[0], lo = pts[0];
    pts.forEach(function (p) { if (p.M_inf > hi.M_inf) hi = p; if (p.M_inf < lo.M_inf) lo = p; });
    return { factor: f, pts: pts,
      M_I_max: max.apply(null, M_I), M_II_max: max.apply(null, M_II),
      M_inf_max: hi.M_inf, x_inf_max: hi.x,
      M_I_min: min.apply(null, M_I), M_II_min: min.apply(null, M_II),
      M_inf_min: lo.M_inf, x_inf_min: lo.x,
      at: function (x) {
        var best = pts[0];
        pts.forEach(function (p) { if (abs(p.x - x) < abs(best.x - x)) best = p; });
        return best;
      },
      // λ(t) 單調 → 極值只在 t₁（M_I）與 ∞（M_inf）；M_II 實際從未達到
      envelope: function (x) {
        var p = this.at(x);
        return [min(p.M_I, p.M_inf), max(p.M_I, p.M_inf)];
      } };
  };
  BC.simpleSpanDLMoment = function (spans, w, x) {
    var acc = 0;
    for (var i = 0; i < spans.length; i++) {
      if (x <= acc + spans[i] + 1e-9) {
        var xp = min(max(x - acc, 0), spans[i]);
        return w * xp * (spans[i] - xp) / 2;
      }
      acc += spans[i];
    }
    return 0;
  };
  BC.spanBySpanDeadLoad = function (spans, w, nPerSpan) {
    nPerSpan = nPerSpan || 20;
    var total = spans.reduce(function (a, b) { return a + b; }, 0),
        n = nPerSpan * spans.length, xs = [], m1 = [], m2 = [], brk = [], acc = 0, i;
    spans.forEach(function (L) { acc += L; brk.push(acc); });
    for (i = 0; i <= n; i++) {
      var x = total * i / n;
      xs.push(x); m1.push(BC.simpleSpanDLMoment(spans, w, x)); m2.push(BC.contDLMoment(spans, w, x));
    }
    // 支承處補點：束制彎矩的折點都在支承，漏掉會讓線性檢查與極值失真
    brk.slice(0, -1).forEach(function (b) {
      if (xs.every(function (v) { return abs(b - v) > 1e-9; })) {
        xs.push(b); m1.push(BC.simpleSpanDLMoment(spans, w, b)); m2.push(BC.contDLMoment(spans, w, b));
      }
    });
    var ord = xs.map(function (v, k) { return k; }).sort(function (a, b) { return xs[a] - xs[b]; });
    return { xs: ord.map(function (k) { return xs[k]; }),
             M_I: ord.map(function (k) { return m1[k]; }),
             M_II: ord.map(function (k) { return m2[k]; }) };
  };
  // 自我驗證：束制彎矩 ΔM 不對應外載 → 支承間必為線性（任一側算錯即掛）
  BC.redistributionIsLinear = function (res, spans, tol) {
    tol = tol == null ? 1e-6 : tol;
    var edges = [0], acc = 0, scale = 0;
    spans.forEach(function (L) { acc += L; edges.push(acc); });
    res.pts.forEach(function (p) { scale = max(scale, abs(p.dM)); });
    scale = scale || 1;
    for (var k = 0; k < edges.length - 1; k++) {
      var a = edges[k], b = edges[k + 1];
      var seg = res.pts.filter(function (p) { return p.x >= a - 1e-9 && p.x <= b + 1e-9; });
      if (seg.length < 3) continue;
      var x0 = seg[0].x, y0 = seg[0].dM, x1 = seg[seg.length - 1].x, y1 = seg[seg.length - 1].dM;
      if (x1 - x0 <= 0) continue;
      for (var j = 0; j < seg.length; j++) {
        var lin = y0 + (y1 - y0) * (seg[j].x - x0) / (x1 - x0);
        if (abs(seg[j].dM - lin) > tol * scale) return false;
      }
    }
    return true;
  };
  // 先簡支張拉後連續：張拉當下 M₂≡0（靜定），連續後由潛變生成 λ·M₂,cont
  BC.prestressM2Redistribution = function (spans, xs, M2cont, dphi, chi, method, phiR) {
    var f = BC.redistributionFactor(dphi, chi, method, phiR),
        inf = M2cont.map(function (m) { return f.lam * m; }), ip = 0;
    for (var i = 1; i < xs.length; i++) if (abs(xs[i] - spans[0]) < abs(xs[ip] - spans[0])) ip = i;
    return { factor: f, xs: xs, M2_cont: M2cont, M2_inf: inf,
             M2_pier_cont: M2cont[ip], M2_pier_inf: inf[ip] };
  };
  // 簡支時張拉的腱：各跨自己的拋物線（連續腱線形只在合龍後才存在，不可拿來算 M₂,cont）
  BC.simpleSpanTendonSegs = function (spans, eMid, eEnd) {
    eEnd = eEnd || 0;
    var segs = [], acc = 0;
    spans.forEach(function (L) {
      segs.push(BC.parabolaSeg(acc, acc + L, acc + L / 2, eMid, -4 * (eMid - eEnd) / (L * L), 'simple'));
      acc += L;
    });
    return segs;
  };
  // 體系轉換套進連續梁包絡：t₁（自重＝簡支）與 ∞（重分配）兩狀態取不利；
  // 恆載對所求彎矩有利時取 γ_min（墩頂正彎矩接頭的關鍵）。AASHTO 5.14.1.4.2 自動滿足。
  BC.stagedEnvelope = function (spans, rows, wDC, lam, o) {
    o = o || {};
    var M2 = o.M2 || null, pss = !!o.psAtSimple,
        gdc = o.gDC == null ? 1.25 : o.gDC, gdw = o.gDW == null ? 1.50 : o.gDW, gll = o.gLL == null ? 1.75 : o.gLL,
        gdcm = o.gDCmin == null ? 0.90 : o.gDCmin, gdwm = o.gDWmin == null ? 0.65 : o.gDWmin;
    return rows.map(function (r, i) {
      var mI = BC.simpleSpanDLMoment(spans, wDC, r.x), mII = r.M_dc, minf = mI + lam * (mII - mI),
          m2 = M2 ? M2[i] : 0, m2t1 = pss ? 0 : m2, m2inf = pss ? lam * m2 : m2;
      var st = {};
      [['t1', mI, m2t1], ['inf', minf, m2inf]].forEach(function (q) {
        var dc = q[1], mm = q[2];
        st[q[0]] = [dc + r.M_dw + r.M_ll_pos + mm, dc + r.M_dw + r.M_ll_neg + mm,
                    BC.factored(dc, r.M_dw, r.M_ll_pos, true, gdc, gdcm, gdw, gdwm, gll) + mm,
                    BC.factored(dc, r.M_dw, r.M_ll_neg, false, gdc, gdcm, gdw, gdwm, gll) + mm];
      });
      return { x: r.x, M_dc_I: mI, M_dc_II: mII, M_dc_inf: minf, M_dw: r.M_dw,
               M_ll_pos: r.M_ll_pos, M_ll_neg: r.M_ll_neg, M2_t1: m2t1, M2_inf: m2inf,
               Ms_pos: max(st.t1[0], st.inf[0]), Ms_neg: min(st.t1[1], st.inf[1]),
               Mu_pos: max(st.t1[2], st.inf[2]), Mu_neg: min(st.t1[3], st.inf[3]),
               gov_pos: st.t1[2] >= st.inf[2] ? 't1' : 'inf', gov_neg: st.t1[3] <= st.inf[3] ? 't1' : 'inf' };
    });
  };
  BC.simpleSpanDLShear = function (spans, w, x, side) {
    var xs = contXs(spans), i = contShearSpan(xs, spans, x, side || 'R');
    return w * (spans[i] / 2 - (x - xs[i]));
  };
  // 束制剪力 λ(V_II − V_I) 為束制彎矩斜率，每跨常數；端支承 t₁ 為簡支 wL/2 > 連續 3wL/8
  BC.stagedShearRow = function (spans, r, wdc, lam) {
    var vI = BC.simpleSpanDLShear(spans, wdc, r.x, r.side), vinf = vI + lam * (r.V_dc - vI);
    var c = [['t1', vI], ['inf', vinf]].map(function (q) {
      return [q[0], q[1], BC.factored(q[1], r.V_dw, r.V_ll_pos, true), BC.factored(q[1], r.V_dw, r.V_ll_neg, false)];
    });
    var p = c[0][2] >= c[1][2] ? c[0] : c[1], n = c[0][3] <= c[1][3] ? c[0] : c[1],
        ctrl = Math.abs(p[2]) >= Math.abs(n[3]) ? p : n;
    return { x: r.x, side: r.side, V_dc: ctrl[1], V_dw: r.V_dw, V_ll_pos: r.V_ll_pos, V_ll_neg: r.V_ll_neg,
             Vu_pos: p[2], Vu_neg: n[3], I: r.I, V_dc_I: vI, V_dc_II: r.V_dc, V_dc_inf: vinf,
             gov_pos: p[0], gov_neg: n[0] };
  };
  // 連續橫隔梁正彎矩接頭（AASHTO 5.14.1.4.9a／.4）：max(Mu⁺, 0.6Mcr)；梁齡 ≥ 90 天簡化 1.2Mcr。
  // Mcr＝f_r·I_g/y_t：總毛斷面、橫隔梁混凝土 f_r 0.63√f'c（5.4.2.6 基本值）、不計預力。
  BC.positiveMomentConnection = function (MuPos, Ig, yt, fcD, ageDays, fr, d, fy, phi) {
    fr = fr == null ? 0.63 * sqrt(fcD) : fr; fy = fy || 420; phi = phi == null ? 0.9 : phi;
    var Mcr = fr * Ig / yt / 1e6, simp = ageDays != null && ageDays >= 90, Mreq, gov;
    if (simp) { Mreq = 1.2 * Mcr; gov = '1.2Mcr'; }
    else if (MuPos > 0.6 * Mcr) { Mreq = MuPos; gov = 'Mu+'; } else { Mreq = 0.6 * Mcr; gov = '0.6Mcr'; }
    return { fr: fr, Mcr: Mcr, Mu_pos: MuPos, simplified: simp, M_req: Mreq, governs: gov,
             As_est: d ? Mreq * 1e6 / (phi * fy * 0.9 * d) : null };
  };
  // 正彎矩接頭之錨定與配置（同 staging；AASHTO 5.14.1.4.9b/c/d，2026-09-22 NLM 核）
  // 5.14.1.4.9c 延伸鋼絞線設計應力：f_psl=(ℓ_dsh−203)/0.840（服務、斷面開裂）、
  //   f_pul=(ℓ_dsh−203)/0.600（強度）。🔴 **式中沒有 d_b**（英制 (ℓ−8)/0.228、/0.163 ksi-in，
  //   換算 6.895/(0.228×25.4)=1/0.840 ✓）；彎折前自梁面外伸 ≥200 mm；debonded 者不得使用。
  // 5.14.1.4.9b 一般鋼筋：錨定依 Art. 5.11，且須伸展**超過支承面內側邊緣**；多支時截斷點
  //   成對錯開且對稱於梁中心線。5.14.1.4.9d：配置對稱、兩側梁鋼筋須能交錯不衝突。
  BC.STRAND_PROJECT_MIN = 200; BC.STRAND_L0 = 203;
  BC.strandStressExtended = function (l) {
    var e = Math.max(0, l - BC.STRAND_L0); return { f_psl: e / 0.840, f_pul: e / 0.600 };
  };
  BC.strandLengthRequired = function (f, limit) {
    return f * (limit === 'service' ? 0.840 : 0.600) + BC.STRAND_L0;
  };
  BC.posConnStrand = function (Mreq, d, lDsh, Astrand, projection, phi, jd) {
    Astrand = Astrand == null ? 98.7 : Astrand; projection = projection == null ? 250 : projection;
    phi = phi == null ? 0.9 : phi; jd = jd == null ? 0.9 : jd;
    var f = BC.strandStressExtended(lDsh), capEach = phi * f.f_pul * Astrand * jd * d / 1e6;
    var nReq = capEach > 0 ? Mreq / capEach : Infinity;
    var nUse = isFinite(nReq) ? Math.ceil(nReq / 2) * 2 : 0;
    return { kind: 'strand', M_req: Mreq, d: d, n_req: nReq, n_use: nUse, A_each: Astrand,
             f_design: f.f_pul, f_psl: f.f_psl, l_dsh: lDsh, phiMn: nUse * capEach,
             ok: nUse * capEach >= Mreq, project_ok: projection >= BC.STRAND_PROJECT_MIN };
  };
  BC.posConnRebar = function (Mreq, d, ld, devAvailable, barArea, fy, phi, jd) {
    barArea = barArea == null ? 387 : barArea; fy = fy || 420;
    phi = phi == null ? 0.9 : phi; jd = jd == null ? 0.9 : jd;
    var capEach = phi * fy * barArea * jd * d / 1e6, nReq = capEach > 0 ? Mreq / capEach : Infinity;
    var nUse = isFinite(nReq) ? Math.ceil(nReq / 2) * 2 : 0;
    return { kind: 'rebar', M_req: Mreq, d: d, n_req: nReq, n_use: nUse, A_each: barArea,
             f_design: fy, ld: ld, dev_available: devAvailable, dev_ok: devAvailable >= ld,
             phiMn: nUse * capEach, ok: nUse * capEach >= Mreq };
  };

  // 多階段逐跨施工（同 staging.multi_stage_redistribution）：各跨材齡不同 → 每階段自有 λ
  //   M(∞) = Σ_i [ M_I,i + λ_i (M_II,i − M_I,i) ]，λ_i 以該階段之 (t0_i, t_c,i) 求嚴格 AAEM。
  // 🔴 便覽 §3.5.4 但書：材齡差 > 365 天時，單一 λ 的簡化不適用（回報 age_gap_warn）。
  // 🔴 簡化的偏差方向**依各階段 M_I／M_II 的正負組合而定**，不可假設單一 λ 偏保守。
  BC.AGE_GAP_LIMIT_DAYS = 365;
  BC.multiStageRedistribution = function (stages, o) {
    o = o || {};
    var chi = o.chi == null ? 0.8 : o.chi, H = o.H == null ? 75 : o.H,
        VS = o.VS == null ? 150 : o.VS, fci = o.fci == null ? 32 : o.fci;
    var rows = [], Mtot = 0, MI = 0, MII = 0, tcs = [];
    stages.forEach(function (st) {
      var sp = BC.stagingPhi(st.t0, st.t_c, H, VS, fci);
      var lam = BC.redistributionFactor(sp.dphi_sub, chi, 'trost', sp.phi_load_t1).lam;
      var Mi = st.M_I + lam * (st.M_II - st.M_I);
      rows.push({ name: st.name, lam: lam, M_I: st.M_I, M_II: st.M_II, M: Mi });
      Mtot += Mi; MI += st.M_I; MII += st.M_II; tcs.push(st.t_c);
    });
    var gap = tcs.length ? Math.max.apply(null, tcs) - Math.min.apply(null, tcs) : 0;
    var t0s = o.single_t0 == null ? stages[0].t0 : o.single_t0,
        tc1 = o.single_tc == null ? Math.max.apply(null, tcs) : o.single_tc;
    var sp1 = BC.stagingPhi(t0s, tc1, H, VS, fci),
        lam1 = BC.redistributionFactor(sp1.dphi_sub, chi, 'trost', sp1.phi_load_t1).lam;
    var den = MII - MI;
    return { rows: rows, M_total: Mtot, M_I_total: MI, M_II_total: MII,
             lam_equiv: Math.abs(den) > 1e-12 ? (Mtot - MI) / den : 0,
             age_gap: gap, age_gap_warn: gap > BC.AGE_GAP_LIMIT_DAYS,
             single_lam: lam1, M_single: MI + lam1 * (MII - MI) };
  };
  // 逐跨施工之各階段彎矩拆解（同 staging.stage_moments_span_by_span）
  // 💡 自我驗證：Σ M_I,k ＝簡支、Σ M_II,k ＝連續（疊加原理）
  BC.stageMomentsSpanBySpan = function (spans, w, x) {
    var xs = contXs(spans), out = [];
    spans.forEach(function (L, k) {
      var m1 = (x >= xs[k] - 1e-9 && x <= xs[k + 1] + 1e-9)
        ? w * (x - xs[k]) * (L - (x - xs[k])) / 2 : 0;
      var wv = spans.map(function (_, j) { return j === k ? w : 0; });
      var M0 = BC.contSupportMomentsUniform(spans, wv), i = contSpanOf(xs, x),
          a = x - xs[i], Li = spans[i];
      out.push([m1, M0[i] + (M0[i + 1] - M0[i]) * a / Li + wv[i] * a * (Li - a) / 2]);
    });
    return out;
  };
  // 逐跨施工排程 → 各跨 (t0, t_c)：第 k 跨於 k·daysPerSpan 天澆置，全橋於最後一跨合龍
  BC.spanBySpanSchedule = function (nSpans, daysPerSpan, t0Offset, firstCastAge) {
    t0Offset = t0Offset == null ? 28 : t0Offset; firstCastAge = firstCastAge == null ? 28 : firstCastAge;
    var last = (nSpans - 1) * daysPerSpan, out = [];
    for (var k = 0; k < nSpans; k++) {
      var t0 = k === 0 ? firstCastAge : t0Offset;
      out.push([t0, Math.max(t0 + (last - k * daysPerSpan), t0)]);
    }
    return out;
  };

  // ── 懸臂工法之體系轉換（C2，同 staging 的 cantilever_*）────────────────────
  // 🔴 符號對照：便覽／H3 卡的 X₁＝M_I（合龍前懸臂體系）、X₀＝M_II（同載重一次完工）
  //    X_final = X₁ + λ(X₀ − X₁)，與 M(∞) = M_I + λ(M_II − M_I) 同式。
  // 🔴 X₀ 與 X₁ 必須是**同一組載重**：把合龍後才施加的載重（支架段落架、SDL）塞進 X₀，
  //    會讓它只拿到 λ 倍而少算 (1−λ) 倍；墩頂負彎矩少算＝**偏不安全**。誤差恰 (λ−1)·M_post。
  BC.cantileverUnits = function (spans, piers, closures) {
    var total = spans.reduce(function (a, b) { return a + b; }, 0), units = [];
    closures.forEach(function (c) { if (c < 0 || c > total) throw new Error('合龍點不在橋長內'); });
    piers.forEach(function (xp) {
      var lf = closures.filter(function (c) { return c < xp - 1e-9; }),
          rt = closures.filter(function (c) { return c > xp + 1e-9; });
      var aL = lf.length ? xp - Math.max.apply(null, lf) : 0,
          aR = rt.length ? Math.min.apply(null, rt) - xp : 0;
      units.push({ x_pier: xp, a_left: aL, a_right: aR, x_tip_left: xp - aL, x_tip_right: xp + aR });
    });
    for (var i = 0; i + 1 < units.length; i++)
      if (units[i].x_tip_right > units[i + 1].x_tip_left + 1e-9) throw new Error('懸臂單元重疊');
    return units;
  };
  // 由跨徑推出各內支承出平衡懸臂之 (piers, closures)：墩間跨合龍於跨中，端跨合龍距墩 endFrac×跨長
  BC.cantileverLayout = function (spans, endFrac) {
    endFrac = endFrac == null ? 0.5 : endFrac;
    var xs = [0]; spans.forEach(function (L) { xs.push(xs[xs.length - 1] + L); });
    var piers = xs.slice(1, -1);
    if (!piers.length) return { piers: [], closures: [] };
    var cl = [piers[0] - endFrac * spans[0]];
    for (var i = 0; i + 1 < piers.length; i++) cl.push((piers[i] + piers[i + 1]) / 2);
    cl.push(piers[piers.length - 1] + endFrac * spans[spans.length - 1]);
    return { piers: piers, closures: cl.sort(function (a, b) { return a - b; }) };
  };
  BC.cantileverCovered = function (units) {
    return units.map(function (u) { return [u.x_tip_left, u.x_tip_right]; });
  };
  BC.cantileverFalsework = function (spans, units) {
    var total = spans.reduce(function (a, b) { return a + b; }, 0), out = [], cur = 0;
    BC.cantileverCovered(units).slice().sort(function (a, b) { return a[0] - b[0]; })
      .forEach(function (iv) {
        if (iv[0] > cur + 1e-9) out.push([cur, iv[0]]);
        cur = max(cur, iv[1]);
      });
    if (cur < total - 1e-9) out.push([cur, total]);
    return out;
  };
  BC.cantileverMI = function (units, w, x, side) {
    side = side || 'R';
    for (var i = 0; i < units.length; i++) {
      var u = units[i];
      if (x >= u.x_tip_left - 1e-9 && x <= u.x_tip_right + 1e-9) {
        if (abs(x - u.x_pier) < 1e-9) {
          var a = side.toUpperCase().charAt(0) === 'R' ? u.a_right : u.a_left;
          return -w * a * a / 2;
        }
        if (x > u.x_pier) return -w * Math.pow(u.x_tip_right - x, 2) / 2;
        return -w * Math.pow(x - u.x_tip_left, 2) / 2;
      }
    }
    return 0;
  };
  // 🔴 overhang：掛籃永遠掛在已澆節塊尖端**之外**，不放寬就會被整個丟掉（初期掛籃項大於節塊自重）
  BC.cantileverMIPoints = function (units, loads, x, side, overhang) {
    side = side || 'R'; overhang = overhang || 0;
    for (var i = 0; i < units.length; i++) {
      var u = units[i];
      if (!(x >= u.x_tip_left - 1e-9 && x <= u.x_tip_right + 1e-9)) continue;
      var right = (x > u.x_pier + 1e-9) ||
                  (abs(x - u.x_pier) < 1e-9 && side.toUpperCase().charAt(0) === 'R');
      var lo = u.x_tip_left - overhang, hi = u.x_tip_right + overhang, m = 0;
      loads.forEach(function (ld) {
        var p = ld[0], P = ld[1];
        if (!(p >= lo - 1e-9 && p <= hi + 1e-9)) return;
        if (right && p >= x - 1e-9) m += -P * (p - x);
        else if (!right && p <= x + 1e-9) m += -P * (x - p);
      });
      return m;
    }
    return 0;
  };
  BC.cantileverUnbalanced = function (units, w, extras, overhang) {
    extras = extras || []; overhang = overhang || 0;
    return units.map(function (u) {
      var mL = BC.cantileverMI([u], w, u.x_pier, 'L') + BC.cantileverMIPoints([u], extras, u.x_pier, 'L', overhang);
      var mR = BC.cantileverMI([u], w, u.x_pier, 'R') + BC.cantileverMIPoints([u], extras, u.x_pier, 'R', overhang);
      return [u.x_pier, mL - mR];
    });
  };
  function cantGrid(spans, units, nPer) {
    var total = spans.reduce(function (a, b) { return a + b; }, 0), n = nPer * spans.length, xs = [], i;
    for (i = 0; i <= n; i++) xs.push(total * i / n);
    var key = [], acc = 0;
    spans.forEach(function (L) { acc += L; key.push(acc); });
    units.forEach(function (u) { key.push(u.x_pier); key.push(u.x_tip_left); key.push(u.x_tip_right); });
    key.forEach(function (k) {
      if (k >= 0 && k <= total && xs.every(function (v) { return abs(k - v) > 1e-9; })) xs.push(k);
    });
    return xs.sort(function (a, b) { return a - b; });
  }
  BC.cantileverDeadLoad = function (spans, units, w, nPer, side) {
    nPer = nPer == null ? 20 : nPer;
    var loads = BC.cantileverCovered(units).map(function (iv) { return [iv[0], iv[1], w]; });
    var xs = cantGrid(spans, units, nPer);
    return { xs: xs,
      M_I: xs.map(function (x) { return BC.cantileverMI(units, w, x, side); }),
      M_II: xs.map(function (x) { return BC.contMomentPartialUDL(spans, loads, x); }) };
  };
  BC.cantileverX1X0 = function (spans, units, w, dphi, o) {
    o = o || {};
    var wSdl = o.w_sdl || 0, nPer = o.n_per_span == null ? 20 : o.n_per_span;
    var dl = BC.cantileverDeadLoad(spans, units, w, nPer, o.side);
    var res = BC.creepRedistribution(dl.xs, dl.M_I, dl.M_II, dphi, o.chi, o.method, o.phi_restraint);
    var lam = res.factor.lam;
    var post = BC.cantileverFalsework(spans, units).map(function (iv) { return [iv[0], iv[1], w]; });
    var mPost = dl.xs.map(function (x) {
      return BC.contMomentPartialUDL(spans, post, x) + (wSdl ? BC.contDLMoment(spans, wSdl, x) : 0);
    });
    var mInf = res.pts.map(function (p) { return p.M_inf; });
    var mTot = mInf.map(function (v, i) { return v + mPost[i]; });
    var mWrong = dl.M_I.map(function (a, i) { return a + lam * (dl.M_II[i] + mPost[i] - a); });
    function pick(x) {
      var bi = 0;
      dl.xs.forEach(function (v, i) { if (abs(v - x) < abs(dl.xs[bi] - x)) bi = i; });
      return [dl.xs[bi], dl.M_I[bi], dl.M_II[bi], mInf[bi]];
    }
    var cls = [];
    units.forEach(function (u) { cls.push(u.x_tip_left); cls.push(u.x_tip_right); });
    cls = cls.filter(function (c, i) { return cls.indexOf(c) === i; })
             .filter(function (c) { return units.every(function (u) { return abs(c - u.x_pier) > 1e-9; }); })
             .sort(function (a, b) { return a - b; });
    var ip = 0;
    dl.xs.forEach(function (v, i) { if (abs(v - units[0].x_pier) < abs(dl.xs[ip] - units[0].x_pier)) ip = i; });
    return { lam: lam, xs: dl.xs, M_I: dl.M_I, M_II: dl.M_II, M_inf: mInf,
             piers: units.map(function (u) { return pick(u.x_pier); }),
             closures: cls.map(pick),
             unbalanced: BC.cantileverUnbalanced(units, w),
             linear_ok: BC.redistributionIsLinear(res, spans),
             M_post: mPost, M_dead_total: mTot, M_lumped_wrong: mWrong,
             lumped_err_pier: mWrong[ip] - mTot[ip] };
  };
  // 節塊施工循環 → 各節塊 (t0, t_c)：第 k 塊（k=0 最靠墩）於 k·days 澆置、合龍日 n·days+extra
  BC.cantileverSchedule = function (nSeg, daysPerSeg, t0, daysToClosure) {
    t0 = t0 == null ? 7 : t0; daysToClosure = daysToClosure || 0;
    var close = nSeg * daysPerSeg + daysToClosure, out = [];
    for (var k = 0; k < nSeg; k++) out.push([t0, max(close - k * daysPerSeg, t0)]);
    return out;
  };
  // 各節塊自重 → StageSpec 串列（餵 multiStageRedistribution）；segments=[[name,p,G]]
  BC.cantileverStages = function (spans, units, segments, x, schedule, side, overhang) {
    var sch = schedule || BC.cantileverSchedule(segments.length, 7);
    if (sch.length !== segments.length) throw new Error('schedule 長度須與 segments 相同');
    return segments.map(function (sg, i) {
      return { name: sg[0], t0: sch[i][0], t_c: sch[i][1],
               M_I: BC.cantileverMIPoints(units, [[sg[1], sg[2]]], x, side, overhang),
               M_II: sg[2] * BC.contMomentIL(spans, x, sg[1]) };
    });
  };

  // 負彎矩接頭（同 staging.negative_moment_connection；AASHTO 5.14.1.4.5/.6/.7/.8＋5.11.1.2.3）
  // 5.14.1.4.5：正、負彎矩接頭**兩者都要設**，不論連續程度。
  // 5.14.1.4.8：橋面板縱向鋼筋依強度極限負彎矩配置；**錨定端須落在強度狀態下受壓之區域**、
  //   截斷點須錯開；橋面板全部縱向筋皆可計入；**無複合橋面板時梁間接頭為必須**（依 5.11.5）。
  // 5.11.1.2.3：**至少 1/3** 負彎矩鋼筋延伸超過反曲點 ≥ max(d, 12d_b, 0.0625×淨跨)。
  // 5.14.1.4.6：梁頂拉應力用 Table 5.9.4.1.2-1，但**以 f'c 代 f'ci**：有握裹 0.63√f'c、
  //   無握裹 min(0.25√f'c, 1.38 MPa)；以 Service III 計算。
  BC.negTopTensionLimit = function (fc, bonded) {
    return bonded === false ? Math.min(0.25 * sqrt(fc), 1.38) : 0.63 * sqrt(fc);
  };
  BC.negativeMomentConnection = function (MuNeg, d, clearSpan, o) {
    o = o || {};
    var barArea = o.bar_area == null ? 387 : o.bar_area, db = o.db == null ? 22.2 : o.db,
        fy = o.fy || 420, phi = o.phi == null ? 0.9 : o.phi, jd = o.jd == null ? 0.9 : o.jd,
        fc = o.fc == null ? 40 : o.fc, bonded = o.bonded !== false,
        composite = o.composite_deck !== false;
    var AsReq = MuNeg * 1e6 / (phi * fy * jd * d), nReq = AsReq / barArea, nUse = Math.ceil(nReq);
    var AsProv = o.As_provided == null ? nUse * barArea : o.As_provided;
    var embedReq = Math.max(d, 12 * db, 0.0625 * clearSpan),
        have = o.embed_beyond_PI == null ? 0 : o.embed_beyond_PI,
        lim = BC.negTopTensionLimit(fc, bonded), st = o.sigma_top;
    return { Mu_neg: MuNeg, d: d, As_req: AsReq, n_req: nReq, n_use: nUse, As_provided: AsProv,
             ok_strength: AsProv >= AsReq - 1e-9, embed_req: embedReq, embed_have: have,
             embed_ok: have >= embedReq - 1e-9, n_one_third: Math.ceil(nUse / 3),
             sigma_top: st == null ? null : st, sigma_limit: lim,
             service_ok: st == null || st <= lim + 1e-9, connection_required: !composite };
  };

  // AASHTO LRFD 5.4.2.3.2 潛變係數（台灣規範無 φ(t) 模式；同 staging.aashto_creep）
  BC.aashtoCreep = function (t, ti, H, VS, fci) {
    H = H == null ? 75 : H; VS = VS == null ? 150 : VS; fci = fci == null ? 32 : fci;
    var ks = max(1, 1.45 - 0.0051 * VS), khc = 1.56 - 0.008 * H, kf = 35 / (7 + fci),
        ktd = t === Infinity ? 1 : t / (61 - 0.58 * fci + t);
    return { ks: ks, khc: khc, kf: kf, ktd: ktd, psi: 1.9 * ks * khc * kf * ktd * Math.pow(ti, -0.118) };
  };
  BC.stagingPhi = function (t0, t1, H, VS, fci) {
    if (t1 < t0) throw new Error('連續材齡不可早於載重材齡');
    var inf0 = BC.aashtoCreep(Infinity, t0, H, VS, fci).psi,
        at1 = t1 > t0 ? BC.aashtoCreep(t1 - t0, t0, H, VS, fci).psi : 0;
    return { t0: t0, t1: t1, phi_inf: inf0, phi_t1: at1, dphi_sub: inf0 - at1,
             phi_load_t1: BC.aashtoCreep(Infinity, t1, H, VS, fci).psi };
  };
  // 箱梁 V/S：C5.4.2.3.2 封閉箱室內周長只計 50%
  BC.boxVolumeSurface = function (topW, topT, botW, botT, webT, nWeb, h, innerFactor) {
    innerFactor = innerFactor == null ? 0.5 : innerFactor;
    var hw = h - topT - botT, A = topW * topT + botW * botT + nWeb * webT * hw,
        pOut = topW + (topW - botW) + 2 * topT + 2 * hw + 2 * botT + botW,
        pIn = 2 * (botW - nWeb * webT) + 2 * max(0, nWeb - 1) * hw;
    return A / (pOut + innerFactor * pIn);
  };
  BC.timingSensitivityAashto = function (t0, t1List, M_I_pier, M_II_pier, H, VS, fci, chi) {
    chi = chi == null ? 0.8 : chi;
    return t1List.map(function (t1) {
      var sp = BC.stagingPhi(t0, t1, H, VS, fci),
          ex = BC.redistributionFactor(sp.dphi_sub, chi, 'trost', sp.phi_load_t1).lam,
          ap = BC.redistributionFactor(sp.dphi_sub, chi, 'trost').lam;
      return { t1: t1, dphi: sp.dphi_sub, phi_r: sp.phi_load_t1, lam_exact: ex, lam_approx: ap,
               M_pier: M_I_pier + ex * (M_II_pier - M_I_pier) };
    });
  };
  BC.timingSensitivity = function (phiInf, rows, M_I_pier, M_II_pier, chi, method) {
    return rows.map(function (r) {
      var d = max(0, phiInf - r[1]), f = BC.redistributionFactor(d, chi, method);
      return { t1_label: r[0], phi_t1: r[1], dphi: d, lam: f.lam,
               M_pier: M_I_pier + f.lam * (M_II_pier - M_I_pier) };
    });
  };

  BC.taiwanLaneReduction = function (lanes) { return lanes <= 2 ? 1.0 : (lanes === 3 ? 0.90 : 0.75); };
  BC.contImpactLength = function (spans, x, sign) {
    var xs = contXs(spans), n = spans.length, i = contSpanOf(xs, x), j;
    if (sign > 0 || n === 1) return spans[i];
    for (j = 1; j < n; j++) if (Math.abs(x - xs[j]) < 1e-6) return (spans[j - 1] + spans[j]) / 2;
    var li = i >= 1, ri = i <= n - 2;
    if (li && (!ri || x - xs[i] <= xs[i + 1] - x)) return (spans[i - 1] + spans[i]) / 2;
    return (spans[i] + spans[i + 1]) / 2;
  };
  function contCache(spans, stiff) {
    var xs = contXs(spans), memo = {}, flex = stiff ? BC.contFlex(spans, stiff.Irel, stiff.breaks) : null;
    var getX = function (p) { var key = Math.round(p * 1e6);
      return memo[key] || (memo[key] = flex ? flex.supportMomentsPoint(p) : BC.contSupportMomentsPoint(spans, p)); };
    return { spans: spans, xs: xs, tot: xs[xs.length - 1], flex: flex, X: getX,
      eta: function (x, p) { return contEta(spans, xs, getX(p), x, p); } };
  }
  // HS20-44 中後軸距 V 可調 4.25～9.15（§3.6）：連續梁掃描 V；簡支恆 4.25 最不利（同 influence_cont.py）
  var TW_GRID = 0.05;
  BC.taiwanRearSpacings = function (step) {
    step = step || 0.25; var out = [], v = 4.25;
    while (v < 9.15 - 1e-9) { out.push(Math.round(v * 1e6) / 1e6); v += step; }
    out.push(9.15); return out;
  };
  function truckSetsV(V) {
    var P = BC.TW.P, d = [0, 4.25, 4.25 + V];
    return [[P, d], [P.slice().reverse(), d.slice().reverse().map(function (t) { return d[2] - t; })]];
  }
  // vals[k]＝η(k·0.05)；atFn(p, side) 另將各軸逐一放在斷面 x 上（左右極限）
  function truckExtremes(vals, tot, x, atFn, spacings, sEvery) {
    sEvery = sEvery || 2; var N = vals.length - 1, bp = 0, bn = 0, vp = 4.25, vn = 4.25;
    (spacings || BC.taiwanRearSpacings()).forEach(function (V) {
      truckSetsV(V).forEach(function (set) {
        var Ps = set[0], di = set[1].map(function (t) { return Math.round(t / TW_GRID); }), sp = di[2], k, j, v, q;
        for (k = -sp; k <= N; k += sEvery) {
          v = 0;
          for (q = 0; q < 3; q++) { j = k + di[q]; if (j >= 0 && j <= N) v += Ps[q] * vals[j]; }
          if (v > bp) { bp = v; vp = V; } if (v < bn) { bn = v; vn = V; }
        }
        if (atFn) for (var a = 0; a < 3; a++) {
          var s0 = x - set[1][a];
          [-1, 1].forEach(function (sd) {
            var w = 0;
            for (var i2 = 0; i2 < 3; i2++) { var ax = s0 + set[1][i2];
              if (ax >= -1e-9 && ax <= tot + 1e-9) w += Ps[i2] * atFn(Math.min(Math.max(ax, 0), tot), i2 === a ? sd : 0); }
            if (w > bp) { bp = w; vp = V; } if (w < bn) { bn = w; vn = V; }
          });
        }
      });
    });
    return [bp, bn, vp, vn];
  }
  function contLiveAt(c, x, step, grid) {
    var spans = c.spans, xs = c.xs, tot = c.tot, A = BC.TW, j, q;
    var N = Math.round(tot / TW_GRID), vals = [];
    for (q = 0; q <= N; q++) vals.push(c.eta(x, Math.min(q * TW_GRID, tot)));
    var te = truckExtremes(vals, tot, x, function (p) { return c.eta(x, p); }), tp = te[0], tn = te[1];
    var posA = 0, negA = 0, spanMin = [], etaMax = 0;
    for (j = 0; j < spans.length; j++) {
      var h = spans[j] / grid, vals = [], mn = Infinity;
      for (q = 0; q <= grid; q++) { var e = c.eta(x, xs[j] + q * h); vals.push(e); if (e < mn) mn = e; if (e > etaMax) etaMax = e; }
      for (q = 0; q < grid; q++) { var u = vals[q], w = vals[q + 1];
        if (u >= 0 && w >= 0) posA += (u + w) * h / 2;
        else if (u <= 0 && w <= 0) negA += (u + w) * h / 2;
        else { var r = u / (u - w);
          if (u > 0) { posA += u * r * h / 2; negA += w * (1 - r) * h / 2; }
          else { negA += u * r * h / 2; posA += w * (1 - r) * h / 2; } } }
      spanMin.push(mn);
    }
    var negs = spanMin.filter(function (v) { return v < 0; }).sort(function (a, b) { return a - b; });
    var lanePos = A.lane * posA + A.PM * etaMax, laneNeg = A.lane * negA + A.PM * ((negs[0] || 0) + (negs[1] || 0));
    return { x: x, truck_pos: tp, truck_neg: tn, lane_pos: lanePos, lane_neg: laneNeg,
             pos: Math.max(tp, lanePos), neg: Math.min(tn, laneNeg),
             I_pos: BC.taiwanImpact(BC.contImpactLength(spans, x, 1)), I_neg: BC.taiwanImpact(BC.contImpactLength(spans, x, -1)),
             V_pos: te[2], V_neg: te[3] };
  }
  BC.taiwanContLiveMoment = function (spans, x, step, grid) { return contLiveAt(contCache(spans), x, step || 0.25, grid || 400); };
  // 回傳各斷面：x、M_dc、M_dw、M_ll_pos/neg（含衝擊×車道數×折減）、Ms/Mu 正負、I_pos/I_neg（不含預力 M2）
  // stiff：變斷面剖面（BC.haunchProfile）→ 變 EI 柔度、自重 × A(x)/A_ref
  BC.taiwanContEnvelope = function (spans, wdc, wdw, lanes, nPerSpan, step, grid, stiff) {
    nPerSpan = nPerSpan || 20; step = step || 0.25; grid = grid || 400;
    var c = contCache(spans, stiff), xs = c.xs, fac = lanes * BC.taiwanLaneReduction(lanes), pts = [0], i, k;
    var wdcFn = stiff ? function (t) { return wdc * stiff.Arel(t); } : null, wdwFn = function () { return wdw; };
    var Xdc = c.flex ? c.flex.supportMomentsDist(wdcFn) : null, Xdw = c.flex ? c.flex.supportMomentsDist(wdwFn) : null;
    for (i = 0; i < spans.length; i++) for (k = 1; k <= nPerSpan; k++) pts.push(xs[i] + spans[i] * k / nPerSpan);
    return pts.map(function (x) {
      var dc = c.flex ? c.flex.momentAt(Xdc, x, c.flex.M0Dist(wdcFn, x)) : BC.contDLMoment(spans, wdc, x);
      var dw = c.flex ? c.flex.momentAt(Xdw, x, c.flex.M0Dist(wdwFn, x)) : BC.contDLMoment(spans, wdw, x);
      var lv = contLiveAt(c, x, step, grid);
      var lp = lv.pos * (1 + lv.I_pos) * fac, ln = lv.neg * (1 + lv.I_neg) * fac;
      return { x: x, M_dc: dc, M_dw: dw, M_ll_pos: lp, M_ll_neg: ln, Ms_pos: dc + dw + lp, Ms_neg: dc + dw + ln,
               Mu_pos: BC.factored(dc, dw, lp, true), Mu_neg: BC.factored(dc, dw, ln, false), I_pos: lv.I_pos, I_neg: lv.I_neg };
    });
  };

  // ── 連續梁剪力影響線與包絡（同 influence_cont.py）── V＝V0＋(X_{i+1}−X_i)/L_i；side 'L'/'R'；車道集中 116 kN；衝擊長度＝至較遠支點
  function contShearSpan(xs, spans, x, side) {
    for (var j = 1; j < spans.length; j++) if (Math.abs(x - xs[j]) < 1e-9) return side === 'L' ? j - 1 : j;
    return contSpanOf(xs, x);
  }
  function contV0Point(xs, spans, i, x, p, ps) {
    var L = spans[i], a = x - xs[i], q = p - xs[i];
    if (q < -1e-9 || q > L + 1e-9) return 0;
    if (Math.abs(q - a) < 1e-9) return ps < 0 ? -q / L : (L - q) / L;
    return q < a ? -q / L : (L - q) / L;
  }
  BC.contShearIL = function (spans, x, p, side, ps, Xp) {
    side = side || 'R'; var xs = contXs(spans), i = contShearSpan(xs, spans, x, side);
    var X = Xp || BC.contSupportMomentsPoint(spans, p);
    var inSpan = p >= xs[i] - 1e-9 && p <= xs[i + 1] + 1e-9;
    var V0 = inSpan ? contV0Point(xs, spans, i, x, p, ps ? ps : (side === 'R' ? 1 : -1)) : 0;
    return V0 + (X[i + 1] - X[i]) / spans[i];
  };
  BC.contDLShear = function (spans, w, x, side, flex, wFn, Xd) {
    side = side || 'R'; var xs = contXs(spans), i = contShearSpan(xs, spans, x, side), L = spans[i], a = x - xs[i], X, V0;
    if (!flex) { X = BC.contSupportMomentsUniform(spans, spans.map(function () { return w; })); V0 = w * (L / 2 - a); }
    else {
      var wf = wFn || function () { return w; };
      X = Xd || flex.supportMomentsDist(wf);
      var RA = BC.gaussIntegrate(function (t) { return wf(t) * (L - (t - xs[i])) / L; }, xs[i], xs[i + 1], flex.breaks || [], 0.25);
      V0 = RA - BC.gaussIntegrate(wf, xs[i], x, flex.breaks || [], 0.25);
    }
    return V0 + (X[i + 1] - X[i]) / L;
  };
  BC.contShearImpactLength = function (spans, x, side) {
    var xs = contXs(spans), i = contShearSpan(xs, spans, x, side || 'R'), a = x - xs[i]; return Math.max(a, spans[i] - a);
  };
  function contLiveShearAt(c, x, side, step, grid) {
    var spans = c.spans, xs = c.xs, tot = c.tot, A = BC.TW, P0 = A.P, D0 = A.x, sp = D0[2], i = contShearSpan(xs, spans, x, side);
    function eta(p, ps) { return BC.contShearIL(spans, x, p, side, ps, c.X(p)); }
    var N = Math.round(tot / TW_GRID), vals = [], j, q;
    for (q = 0; q <= N; q++) { var pq = Math.min(q * TW_GRID, tot); vals.push(eta(pq, Math.abs(pq - x) < 1e-9 ? (pq > x ? 1 : -1) : 0)); }
    var te = truckExtremes(vals, tot, x, function (p, sd) { return eta(p, sd); }), tp = te[0], tn = te[1];
    var posA = 0, negA = 0, eMax = 0, eMin = 0;
    for (j = 0; j < spans.length; j++) {
      var a0 = xs[j], b0 = xs[j + 1], L = spans[j], pieces = [[a0, b0]];
      if (j === i && x > a0 && x < b0) pieces = [[a0, x], [x, b0]];
      pieces.forEach(function (pc) {
        var u0 = pc[0], v0 = pc[1], n = Math.max(1, Math.round(grid * (v0 - u0) / L)), h = (v0 - u0) / n, vals = [];
        for (q = 0; q <= n; q++) { var p = u0 + q * h, at = j === i && Math.abs(p - x) < 1e-9;
          vals.push(eta(p, at ? (q === n ? -1 : 1) : 0)); }
        for (q = 0; q < n; q++) { var u = vals[q], w = vals[q + 1];
          if (u >= 0 && w >= 0) posA += (u + w) * h / 2;
          else if (u <= 0 && w <= 0) negA += (u + w) * h / 2;
          else { var r = u / (u - w);
            if (u > 0) { posA += u * r * h / 2; negA += w * (1 - r) * h / 2; }
            else { negA += u * r * h / 2; posA += w * (1 - r) * h / 2; } } }
        vals.forEach(function (v) { if (v > eMax) eMax = v; if (v < eMin) eMin = v; });
      });
    }
    var lp = A.lane * posA + A.PV * eMax, ln = A.lane * negA + A.PV * eMin;
    return { x: x, side: side, truck_pos: tp, truck_neg: tn, lane_pos: lp, lane_neg: ln,
             pos: Math.max(tp, lp), neg: Math.min(tn, ln), I: BC.taiwanImpact(BC.contShearImpactLength(spans, x, side)) };
  }
  BC.taiwanContLiveShear = function (spans, x, side, step, grid, stiff) {
    return contLiveShearAt(contCache(spans, stiff), x, side || 'R', step || 0.25, grid || 400);
  };
  BC.taiwanContShearAt = function (spans, x, side, wdc, wdw, lanes, step, grid, stiff, cache) {
    var c = cache || contCache(spans, stiff), fac = lanes * BC.taiwanLaneReduction(lanes), dc, dw;
    if (!c.flex) { dc = BC.contDLShear(spans, wdc, x, side); dw = BC.contDLShear(spans, wdw, x, side); }
    else {
      var wf = function (t) { return wdc * stiff.Arel(t); };
      dc = BC.contDLShear(spans, wdc, x, side, c.flex, wf); dw = BC.contDLShear(spans, wdw, x, side, c.flex);
    }
    var lv = contLiveShearAt(c, x, side, step || 0.25, grid || 400), lp = lv.pos * (1 + lv.I) * fac, ln = lv.neg * (1 + lv.I) * fac;
    return { x: x, side: side, V_dc: dc, V_dw: dw, V_ll_pos: lp, V_ll_neg: ln,
             Vu_pos: BC.factored(dc, dw, lp, true), Vu_neg: BC.factored(dc, dw, ln, false), I: lv.I };
  };
  BC.taiwanContShearEnvelope = function (spans, points, wdc, wdw, lanes, stiff) {
    var c = contCache(spans, stiff);
    return points.map(function (pt) { return BC.taiwanContShearAt(spans, pt[0], pt[1], wdc, wdw, lanes, 0.25, 400, stiff, c); });
  };
  // 台灣 §8.20.3 箍筋最大間距（同 shear.stirrup_max_spacing_TW）→ [s_max, 減半, 斷面足夠]
  BC.STD_STIRRUP_SPACINGS = [100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 600];
  BC.stirrupMaxSpacingTW = function (Vs, fc, bw, dv, h) {
    var lim1 = 0.33 * Math.sqrt(fc) * bw * dv, lim2 = 0.66 * Math.sqrt(fc) * bw * dv, halved = Vs > lim1;
    return [halved ? Math.min(0.375 * h, 300) : Math.min(0.75 * h, 600), halved, Vs <= lim2];
  };
  BC.stirrupPickSpacing = function (sAllow, list) {
    var ok = (list || BC.STD_STIRRUP_SPACINGS).filter(function (v) { return v <= sAllow + 1e-9; });
    return ok.length ? Math.max.apply(null, ok) : null;
  };
  BC.stirrupZones = function (xs, picks, supports) {
    var segs = [], sup = supports || [], i;
    for (i = 0; i < xs.length - 1; i++) {
      if (sup.some(function (sx) { return xs[i] < sx && sx < xs[i + 1]; })) continue;
      var a = picks[i], b = picks[i + 1];
      segs.push([xs[i], xs[i + 1], (a == null || b == null) ? null : Math.min(a, b)]);
    }
    sup.forEach(function (sx) {
      var L = -1, Rr = -1;
      xs.forEach(function (x, j) { if (x < sx) L = j; if (x > sx && Rr < 0) Rr = j; });
      if (Rr >= 0 && (L < 0 || xs[L] < sx)) segs.push([sx, xs[Rr], picks[Rr]]);
      if (L >= 0 && (Rr < 0 || xs[Rr] > sx)) segs.push([xs[L], sx, picks[L]]);
    });
    segs.sort(function (p, q) { return p[0] - q[0] || p[1] - q[1]; });
    var out = [];
    segs.forEach(function (sg) {
      if (sg[1] - sg[0] < 1e-9) return;
      var last = out[out.length - 1];
      if (last && Math.abs(last[1] - sg[0]) < 1e-9 && last[2] === sg[2]) last[1] = sg[1]; else out.push(sg.slice());
    });
    return out;
  };
  // TendonGroup 陣列 → prestress_at(x)＝[P kN, e_eff mm, Σ P·de_loc/dx kN]（同 continuous.groups_prestress_at）
  BC.groupsPrestressAt = function (groups, stiff) {
    return function (x) {
      var P = 0, Pe = 0, Ps = 0;
      groups.forEach(function (g) {
        for (var k = 0; k < g.segs.length; k++) { var sg = g.segs[k];
          if (x >= sg.x1 - 1e-9 && x <= sg.x2 + 1e-9) { var p = typeof g.P === 'function' ? g.P(x) : g.P;
            P += p; Pe += p * BC.segE(sg, x); Ps += p * 2 * sg.c * (x - sg.xv) / 1000; return; } }
      });
      if (stiff && P) Ps -= P * (stiff.dyb(x + 0.01) - stiff.dyb(x - 0.01)) / 0.02 / 1000;
      return [P, P ? Pe / P : 0, Ps];
    };
  };
  // 連續梁全長剪力設計掃描＋箍筋分區（同 influence_cont.cont_shear_design_scan）
  BC.contShearDesignScan = function (spans, prestressAt, fm, wdc, wdw, lanes, h, yb, fc, bw, nWebs, Av, step, extra, stiff, sec, fsy, phi, stageLam) {
    step = step || 1; extra = extra || []; fsy = fsy || 420; phi = phi || 0.85;
    var xs = contXs(spans), ref = sec || BC.section(1, 1, yb, h);
    function dvAt(xf) {
      var dv = 0.72 * h;
      for (var it = 0; it < 3; it++) {
        var xx = typeof xf === 'function' ? xf(dv) : xf, q = prestressAt(xx), yc = yb - q[1];
        var neg = BC.contDLMoment(spans, 1, xx) < 0, dp = neg ? yc : h - yc;
        dv = Math.max(0.9 * dp, 0.72 * h);
      }
      return dv;
    }
    var pts = [];
    spans.forEach(function (L, k) {
      var a = xs[k], b = xs[k + 1];
      var dva = dvAt(function (d) { return a + d / 1000; }), dvb = dvAt(function (d) { return b - d / 1000; });
      var x0 = a + dva / 1000, x1 = b - dvb / 1000, n = Math.max(1, Math.round((x1 - x0) / step)), cand = [], seen = {};
      for (var i = 0; i <= n; i++) cand.push(x0 + (x1 - x0) * i / n);
      extra.forEach(function (e) { if (e > x0 && e < x1) cand.push(e); });
      cand.map(function (v) { return Math.round(v * 1e6) / 1e6; }).sort(function (p, q) { return p - q; })
        .forEach(function (x) { if (!seen[x]) { seen[x] = 1; pts.push([x, x - a < b - x ? 'R' : 'L']); } });
    });
    var env = BC.taiwanContShearEnvelope(spans, pts, wdc, wdw, lanes, stiff), AvsMin = BC.AvSminTW(fc, bw, fsy);
    // 逐跨施工：恆載剪力改 t₁（簡支）／∞（重分配）兩狀態取不利（同 staging.staged_shear_row）
    if (stageLam != null) env = env.map(function (r) { return BC.stagedShearRow(spans, r, wdc, stageLam); });
    var rows = pts.map(function (pt, i) {
      var x = pt[0], side = pt[1], r = env[i], dv = dvAt(x), q = prestressAt(x), P = q[0];
      var V2 = BC.secondaryShear(fm, spans, x, side), Vu = BC.designShearWithV2(r.Vu_pos, r.Vu_neg, V2);
      var secx = stiff ? stiff.sectionAt(x) : ref;
      var sh = BC.shearWebAt(P * 1e3, secx, P ? q[2] / P : 0, fc, bw, dv, Vu / nWebs * 1e3, nWebs, phi, fsy);
      var mc = BC.stirrupMaxSpacingTW(Math.max(sh.Vs_req, 0), fc, bw, dv, h);
      var sStr = sh.Av_s_req > 0 ? Av / sh.Av_s_req : Infinity, sAllow = Math.min(sStr, Av / AvsMin, mc[0]);
      return { x: x, side: side, dv: dv, row: r, V2: V2, Vu: Vu, Vu_web: Math.abs(Vu) / nWebs, P: P, Vp_web: sh.Vp / 1e3,
               sigma1: sh.sigma1, Vcw_web: sh.Vcw / 1e3, Av_s_req: sh.Av_s_req, s_code_max: mc[0], halved: mc[1],
               crush_ok: mc[2], s_allow: sAllow, s_pick: mc[2] ? BC.stirrupPickSpacing(sAllow) : null };
    });
    return { rows: rows, zones: BC.stirrupZones(rows.map(function (r) { return r.x; }), rows.map(function (r) { return r.s_pick; }), xs) };
  };
  BC.secondaryShear = function (fm, spans, x, side) {
    var xs = contXs(spans), i = contShearSpan(xs, spans, x, side || 'R'); return (fm.X[i + 1] - fm.X[i]) / spans[i];
  };
  // Vu＝max(|V_載重|, |V_載重+V2|)（次剪力有利時不折減）
  BC.designShearWithV2 = function (VuPos, VuNeg, V2) {
    return [VuPos, VuNeg, VuPos + V2, VuNeg + V2].reduce(function (a, b) { return Math.abs(b) > Math.abs(a) ? b : a; });
  };

  // ── 連續梁鋼腱線形（分段拋物線）＋ 次彎矩 M2（力法）──────────────
  // 慣例：x m、e mm（形心下為正）、P kN、M kN·m（正彎矩為正）；M1 = −P·e/1000
  BC.parabolaSeg = function (x1, x2, xv, ev, c, kind) {
    return { x1: x1, x2: x2, xv: xv, ev: ev, c: c, kind: kind || '',
             R: Math.abs(c) > 1e-12 ? 1e6 / (2 * Math.abs(c)) : Infinity };
  };
  BC.segE = function (g, x) { var d = x - g.xv; return g.ev + g.c * d * d; };
  // G1 §五之補・FHWA Eqn. 3.30–3.35：墩頂段頂點在墩心（長 kPier×跨長）、跨中段頂點在低點，坡度連續
  BC.contTendonSegs = function (spans, eEnd, eMid, ePier, kPier, rEnd) {
    kPier = kPier == null ? 0.12 : kPier; rEnd = rEnd == null ? 0.42 : rEnd;
    var n = spans.length, xs = [0], i;
    spans.forEach(function (L) { xs.push(xs[xs.length - 1] + L); });
    var lows = spans.map(function (L, i) {
      if (n === 1) return xs[0] + L / 2;
      if (i === 0) return xs[0] + rEnd * L;
      if (i === n - 1) return xs[n - 1] + (1 - rEnd) * L;
      return xs[i] + L / 2;
    });
    var segs = [], infl = [];
    function push(x1, x2, xv, ev, c, kind) { if (x2 > x1 + 1e-9) segs.push(BC.parabolaSeg(x1, x2, xv, ev, c, kind)); }
    function halfToPier(xl, xp, L, l2r) {
      var s = Math.abs(xp - xl), b1 = Math.min(kPier * L, s * 0.9), b2 = s - b1, D = eMid - ePier;
      var h1 = D * b1 / s, h2 = D * b2 / s;
      if (l2r) { push(xl, xp - b1, xl, eMid, -h2 / (b2 * b2), 'sag'); push(xp - b1, xp, xp, ePier, h1 / (b1 * b1), 'hog'); infl.push(xp - b1); }
      else { push(xp, xp + b1, xp, ePier, h1 / (b1 * b1), 'hog'); push(xp + b1, xl, xl, eMid, -h2 / (b2 * b2), 'sag'); infl.push(xp + b1); }
    }
    for (i = 0; i < n; i++) {
      var L = spans[i], xa = xs[i], xb = xs[i + 1], xl = lows[i], be, he;
      if (i === 0) { be = xl - xa; he = eMid - eEnd; push(xa, xl, xl, eMid, -he / (be * be), 'sag'); }
      else halfToPier(xl, xa, L, false);
      if (i === n - 1) { be = xb - xl; he = eMid - eEnd; push(xl, xb, xl, eMid, -he / (be * be), 'sag'); }
      else halfToPier(xl, xb, L, true);
    }
    function eAt(x) {
      for (var k = 0; k < segs.length; k++) { var g = segs[k]; if (x >= g.x1 - 1e-9 && x <= g.x2 + 1e-9) return BC.segE(g, x); }
      return BC.segE(segs[segs.length - 1], x);
    }
    return { segs: segs, eAt: eAt, sx: xs, lows: lows, infl: infl, total: xs[n], kPier: kPier, rEnd: rEnd,
             Rmin: Math.min.apply(null, segs.map(function (g) { return g.R; })) };
  };
  // groups: [{P: kN 或 function(x), segs:[parabolaSeg]}]；probe＝判定所在段的里程（錨碇跳躍取單側極限）
  BC.primaryMomentAt = function (groups, x, probe) {
    var t = probe == null ? x : probe, M = 0;
    groups.forEach(function (g) {
      for (var k = 0; k < g.segs.length; k++) { var s = g.segs[k];
        if (t >= s.x1 - 1e-9 && t <= s.x2 + 1e-9) {
          var P = typeof g.P === 'function' ? g.P(x) : g.P;
          M += -P * BC.segE(s, x) / 1000; return; } }
    });
    return M;
  };
  BC.groupBreaks = function (groups) {
    var o = [];
    groups.forEach(function (g) { g.segs.forEach(function (s) { o.push(s.x1, s.x2); }); });
    return o;
  };
  // 力法（EI 常數）：Σ_j F_ij X_j = −∫M1·m_i dx，F_ii=(L左+L右)/3、F_i,i+1=L/6。
  // 子區間對齊折點＋Simpson：M1 分段二次 × m_i 分段線性＝分段三次 → 精確。
  BC.secondaryMomentsForce = function (spans, M1, breaks, nSub, Irel) {
    nSub = nSub || 8; breaks = breaks || [];
    var n = spans.length, sx = [0], i, j, k;
    spans.forEach(function (L) { sx.push(sx[sx.length - 1] + L); });
    var X = sx.map(function () { return 0; });
    function M2at(x) {
      for (var q = 0; q < sx.length - 1; q++) if (x >= sx[q] - 1e-9 && x <= sx[q + 1] + 1e-9) {
        var t = (x - sx[q]) / (sx[q + 1] - sx[q]); return X[q] * (1 - t) + X[q + 1] * t; }
      return 0;
    }
    if (n < 2) return { sx: sx, X: X, F: [], b: [], M2At: M2at };
    var r9 = function (v) { return Math.round(v * 1e9) / 1e9; }, seen = {}, pts = [];
    sx.concat(breaks.filter(function (v) { return v > sx[0] && v < sx[n]; })).forEach(function (v) {
      v = r9(v); if (!seen[v]) { seen[v] = 1; pts.push(v); } });
    pts.sort(function (a, b) { return a - b; });
    function mHat(i, x) {
      var a = sx[i - 1], c0 = sx[i], b = sx[i + 1];
      if (x >= a && x <= c0) return (x - a) / (c0 - a);
      if (x >= c0 && x <= b) return (b - x) / (b - c0);
      return 0;
    }
    var nu = n - 1, bv = [];
    for (i = 0; i < nu; i++) bv.push(0);
    for (k = 0; k < pts.length - 1; k++) {
      var a = pts[k], b = pts[k + 1]; if (b - a < 1e-9) continue;
      var mid = 0.5 * (a + b), h = (b - a) / (2 * nSub);
      for (i = 0; i < nu; i++) {
        if (!(mid >= sx[i] - 1e-9 && mid <= sx[i + 2] + 1e-9)) continue;
        var s = 0;
        for (j = 0; j <= 2 * nSub; j++) {
          var x = a + j * h, w = (j === 0 || j === 2 * nSub) ? 1 : (j % 2 ? 4 : 2);
          s += w * M1(x, mid) * mHat(i + 1, x) / (Irel ? Irel(x) : 1);
        }
        bv[i] += s * h / 3;
      }
    }
    var F = [];
    if (Irel) F = BC.contFlex(spans, Irel, breaks, Math.max.apply(null, spans)).F;   // 變斷面：數值柔度
    else {
      for (i = 0; i < nu; i++) { F.push([]); for (j = 0; j < nu; j++) F[i].push(0); }
      for (i = 0; i < nu; i++) { F[i][i] = (spans[i] + spans[i + 1]) / 3; if (i + 1 < nu) F[i][i + 1] = F[i + 1][i] = spans[i + 1] / 6; }
    }
    var A = F.map(function (row, i) { return row.slice().concat([-bv[i]]); });
    for (i = 0; i < nu; i++) for (j = i + 1; j < nu; j++) {
      var f = A[j][i] / A[i][i]; for (k = i; k <= nu; k++) A[j][k] -= f * A[i][k]; }
    var Xi = []; for (i = 0; i < nu; i++) Xi.push(0);
    for (i = nu - 1; i >= 0; i--) { var sm = A[i][nu]; for (k = i + 1; k < nu; k++) sm -= A[i][k] * Xi[k]; Xi[i] = sm / A[i][i]; }
    for (i = 0; i < nu; i++) X[i + 1] = Xi[i];
    return { sx: sx, X: X, F: F, b: bv, M2At: M2at };
  };
  // 雙系統配束的墩頂局部腱（頂板腱）：每一內支承兩側各一段拋物線，頂點在墩心；b＝min(length, 0.5×該側跨長)
  BC.pierCapTendonSegs = function (spans, length, eAnchor, ePier) {
    var xs = [0], segs = [], lens = [], anchors = [];
    spans.forEach(function (L) { xs.push(xs[xs.length - 1] + L); });
    for (var j = 1; j < spans.length; j++) {
      var xp = xs[j], bl = Math.min(length, 0.5 * spans[j - 1]), br = Math.min(length, 0.5 * spans[j]);
      segs.push(BC.parabolaSeg(xp - bl, xp, xp, ePier, (eAnchor - ePier) / (bl * bl), 'hog'));
      segs.push(BC.parabolaSeg(xp, xp + br, xp, ePier, (eAnchor - ePier) / (br * br), 'hog'));
      lens.push(Math.min(bl, br)); anchors.push(xp - bl, xp + br);
    }
    return { segs: segs, lengths: lens, anchors: anchors,
             Rmin: segs.length ? Math.min.apply(null, segs.map(function (g) { return g.R; })) : Infinity };
  };
  // 錨具滑移損失（公式卡_後張預力短期損失 §二簡化法；同 tendon_profile.anchor_slip_loss）
  // L_set＝√(Δa·Ep/(1000p))；超過 L_m（單端＝全長、雙端＝半長）時損失圖整體平移使 ∫Δσ＝Δa·Ep/1000
  BC.anchorSlipLoss = function (d, Lm, p, slip, Ep) {
    slip = slip == null ? 6 : slip; Ep = Ep || 195000;
    var A = slip * Ep / 1000;
    if (p <= 0) { var c0 = Lm > 0 ? A / Lm : 0; return { dsigma: c0, L_set: Lm, dsigma_anchor: c0, p: p, capped: true }; }
    var Ls = Math.sqrt(A / p);
    if (Ls <= Lm) return { dsigma: d < Ls ? 2 * p * (Ls - d) : 0, L_set: Ls, dsigma_anchor: 2 * p * Ls, p: p, capped: false };
    var c = (A - p * Lm * Lm) / Lm;
    return { dsigma: 2 * p * (Lm - Math.min(d, Lm)) + c, L_set: Lm, dsigma_anchor: 2 * p * Lm + c, p: p, capped: true };
  };
  // 墩頂局部腱（雙端張拉）逐點有效預力：摩擦（自身線形）＋滑移＋其他損失（同 tendon_profile.pier_cap_tendon_force）
  BC.pierCapTendonForce = function (pc, x, fpj, Aps, other, mu, K, slip, Ep) {
    for (var k = 0; k < pc.segs.length; k += 2) {
      var sl = pc.segs[k], sr = pc.segs[k + 1], xa = sl.x1, xb = sr.x2;
      if (x >= xa - 1e-9 && x <= xb + 1e-9) {
        var Lt = xb - xa, local = [[0, sl.x2 - xa, 2 * Math.abs(sl.c) / 1000], [sl.x2 - xa, Lt, 2 * Math.abs(sr.c) / 1000]];
        var s = Math.min(Math.max(x - xa, 0), Lt), r = BC.frictionAt(s, Lt, local, mu, K, 'both');
        var Lm = Lt / 2, rMid = BC.frictionAt(Lm, Lt, local, mu, K, s <= Lm ? 'start' : 'end'), p = Lm > 0 ? fpj * rMid / Lm : 0;
        var sLoss = BC.anchorSlipLoss(Math.min(s, Lt - s), Lm, p, slip, Ep).dsigma, fr = fpj * r, fpe = fpj - fr - sLoss - other;
        return { P: fpe * Aps / 1000, fpe: fpe, friction: fr, slip: sLoss, other: other, s: s, L_t: Lt };
      }
    }
    return { P: 0, fpe: 0, friction: 0, slip: 0, other: other, s: 0, L_t: 0 };
  };
  // 頂板腱墩頂斷面構造檢核（同 tendon_profile.top_slab_tendon_check）
  BC.topSlabTendonCheck = function (h, yb, topT, ePier, n, width, od, cover, dAgg, rule) {
    od = od == null ? 100 : od; cover = cover == null ? 40 : cover;
    var y = yb - ePier, eLo = yb - (h - cover - od / 2), eHi = yb - (h - topT + cover + od / 2);
    var topCov = h - (y + od / 2), botCov = (y - od / 2) - (h - topT);
    var inSlab = topCov >= cover - 1e-6 && botCov >= cover - 1e-6;
    n = Math.max(1, n | 0);
    var pitch = width / n, xo = [], i;
    for (i = 0; i < n; i++) xo.push(-width / 2 + pitch * (i + 0.5));
    var sClear = n > 1 ? pitch - od : Infinity, sReq = BC.ductSpacingRequired(od, dAgg, rule), sOk = sClear >= sReq - 1e-6;
    return { e_hi: eHi, e_lo: eLo, y_center: y, in_slab: inSlab, top_cover: topCov, bot_cover: botCov,
             x_offsets: xo, s_clear: sClear, s_req: sReq, s_ok: sOk, ok: inSlab && sOk };
  };
  // stiff：變斷面剖面（BC.haunchProfile）→ M1＝Σ −P·(e−Δȳ)/1000、被積函數 ÷I_rel
  BC.continuousPrestress = function (spans, groups, nSub, stiff) {
    if (!stiff) return BC.secondaryMomentsForce(spans, function (x, probe) { return BC.primaryMomentAt(groups, x, probe); },
                                                BC.groupBreaks(groups), nSub);
    return BC.secondaryMomentsForce(spans, function (x, probe) {
      var M = 0, t = probe == null ? x : probe;
      groups.forEach(function (g) {
        for (var k = 0; k < g.segs.length; k++) { var sg = g.segs[k];
          if (t >= sg.x1 - 1e-9 && t <= sg.x2 + 1e-9) { var P = typeof g.P === 'function' ? g.P(x) : g.P;
            M += -P * (BC.segE(sg, x) - stiff.dyb(x)) / 1000; return; } }
      });
      return M;
    }, BC.groupBreaks(groups).concat(stiff.breaks), nSub, stiff.Irel);
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
