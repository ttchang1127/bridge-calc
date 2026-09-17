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
  //   'od'         再取 max(…, 孔道外徑)——公式卡 B3 §4-1 台灣欄另列此條，與 AASHTO
  //                「孔道 OD > 100 mm 時淨距 ≥ 孔道外徑」同義；兩卡不一致、原條文未查證，
  //                 故做成可選，'od' 為保守側。
  // φ100、腹板 350、保護層 40 時：'tw' → 兩列（需求 40）、'od' → 單列（需求 100）。
  BC.ductSpacingRequired = function (od, dAgg, rule) {
    dAgg = dAgg == null ? 25 : dAgg;
    var b = Math.max(40, 1.5 * dAgg);
    return rule === 'od' ? Math.max(b, od) : b;
  };
  // 拋物線最小曲率半徑 R = L²/(8a)（mm）；台灣金屬波形管 ≥6,000、日本 ≥100×管道外徑
  BC.radiusOfCurvature = function (a, L) { return L * L / (8 * a); };
  BC.ductLayout = function (nTendons, nWeb, yCgs, o) {
    o = o || {};
    var od = o.od == null ? 100 : o.od, webT = o.webT == null ? 350 : o.webT,
        cover = o.cover == null ? 40 : o.cover, sv = o.sv == null ? 40 : o.sv,
        dAgg = o.dAgg == null ? 25 : o.dAgg, yb = o.yb, h = o.h, rule = o.rule || 'tw';
    nWeb = Math.max(1, nWeb | 0); var n = Math.max(1, nTendons | 0);
    var base = Math.floor(n / nWeb), rem = n % nWeb, perWeb = [], i, j;
    for (i = 0; i < nWeb; i++) perWeb.push(base + (i < rem ? 1 : 0));
    var nPerWeb = Math.max.apply(null, perWeb);

    var sReq = BC.ductSpacingRequired(od, dAgg, rule), avail = webT - 2 * cover;
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
             eMaxPoint: yb == null ? null : yb - cover - od / 2 };
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
  BC.contMomentIL = function (spans, x, p) { var xs = contXs(spans); return contEta(spans, xs, BC.contSupportMomentsPoint(spans, p), x, p); };
  BC.contDLMoment = function (spans, w, x) {
    var xs = contXs(spans), X = BC.contSupportMomentsUniform(spans, spans.map(function () { return w; }));
    var i = contSpanOf(xs, x), a = x - xs[i];
    return contMomentAt(spans, xs, X, x, w * a * (spans[i] - a) / 2);
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
               Mu_pos: 1.25 * dc + 1.5 * dw + 1.75 * lp, Mu_neg: 1.25 * dc + 1.5 * dw + 1.75 * ln, I_pos: lv.I_pos, I_neg: lv.I_neg };
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
             Vu_pos: 1.25 * dc + 1.5 * dw + 1.75 * lp, Vu_neg: 1.25 * dc + 1.5 * dw + 1.75 * ln, I: lv.I };
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
  BC.contShearDesignScan = function (spans, prestressAt, fm, wdc, wdw, lanes, h, yb, fc, bw, nWebs, Av, step, extra, stiff, sec, fsy, phi) {
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
