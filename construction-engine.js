/* construction-engine.js — construction.py + launching.py + segmental.py 的 JS 移植。
   施工階段：H1/H2 支架施拉應力歷程・H3 平衡懸臂・H4 推進 ILM・H5/H6 預鑄節塊。
   對 golden construction_stage_H1H2 / cantilever_H3 / launching_H4 / segmental_H5H6 驗證。
   H1/H2 複用 window.BC.stresses（與箱梁主引擎同符號慣例，壓為負）。
   單位：力 N/kN（標於介面）、彎矩 kN·m、應力 MPa、長度 m 或 mm（標於介面）。 */
(function (global) {
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
    var s = global.BC.stresses(P, sec, e, MswKNm);
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

  if (typeof module !== 'undefined' && module.exports) module.exports = CE;
  global.CE = CE;
})(typeof window !== 'undefined' ? window : globalThis);
