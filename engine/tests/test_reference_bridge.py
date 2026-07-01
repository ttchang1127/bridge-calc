"""回歸測試：40m 參考橋（台灣 HS20-44、2 設計車道、8組×19股最小設計）黃金答案。

黃金答案來源（知識庫，2 車道修正後自洽）：
  算例_40m參考橋活載基準統一 / 算例_後張箱梁服務性應力驗算 / 算例_40m參考橋載重組合。
執行：`pytest`，或 `python test_reference_bridge.py`（無 pytest 亦可，會印重現摘要）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bridgecalc import (Section, Tendon, compute_losses, combinations,
                        lane_live_load, stresses, Pe_min_zero_tension, allowables,
                        shear_web, phiVn, Av_s_min_TW, flexural_strength,
                        deflection_analysis, il_moment_peak, il_shear_simple,
                        abs_max_moment, lane_moment_simple, hl93_per_lane_moment,
                        moment_envelope_simple, taiwan_per_lane_moment, taiwan_per_lane_shear,
                        taiwan_impact, taiwan_truck_moment, taiwan_lane_moment,
                        fatigue_check, stirrup_fatigue,
                        torsion_check, slab_flexure, As_min_slab, temp_gradient_AASHTO,
                        bearing_check, anchorage_check, spiral_local_bearing, expansion_joint,
                        ThermalBand, self_equilibrating_stress, thermal_service_check,
                        secondary_moment, primary_moment, flexural_strength_T, pier_service_stress,
                        tendon_profile, general_zone_burst, node_capacity, f_cu,
                        grout_qc_check, rebar_stress_limit, pc_fatigue_limit,
                        design_life, GROUT, batched_transfer, stage_stress,
                        transfer_tension_limit, variable_depth, cantilever_moment,
                        long_term_deflection, launching_cantilever_moment,
                        launching_span_moment, centric_prestress_required,
                        launching_bottom_stress, n_tendons, jacking_force, bearing_stress,
                        segment_weight, joint_min_prestress, joint_compression,
                        shear_key_design_capacity, shear_key_utilization, bonded_pt_ratio,
                        min_tendon_groups, required_drape, min_section_modulus_Sb)

# ── 40m 參考橋輸入 ──
sec = Section(A=5.065e6, I=3.287e12, yb=1329, h=2100)
ten = Tendon(n_tendons=8, strands_per=19, e=1109)   # 台灣 HS20-44 最小設計
M_DC, M_DW = 24800, 4000
M_LL_IM = lane_live_load(taiwan_per_lane_moment(40), 2, 1.0)   # HS20-44 2 車道 → 6,837


def _close(a, b, tol):
    assert abs(a - b) <= tol, f"{a:.2f} != {b:.2f} (±{tol})"


def test_live_load():
    _close(M_LL_IM, 6837, 5)


def test_losses_nonlinear_coupling():
    L = compute_losses(ten, sec, M_DC, M_DW)
    _close(L.fcgp, 8.60, 0.05)
    _close(L.ES, 28, 1)
    _close(L.creep, 94, 2)
    _close(L.loss_pct * 100, 20.1, 0.3)
    _close(L.fpe, 1115, 3)
    _close(L.Pe / 1e3, 23720, 50)        # kN


def test_combinations():
    c = combinations(M_DC, M_DW, M_LL_IM)
    _close(c["Strength_I"], 48965, 30)
    _close(c["Service_I"], 35637, 5)
    _close(c["Service_III"], 34270, 5)


def test_service_stress():
    L = compute_losses(ten, sec, M_DC, M_DW)
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    st, sb = stresses(L.Pe, sec, ten.e, M_serv)
    _close(sb, -0.91, 0.06)              # 底緣全壓（HS20 較輕，餘裕大）
    _close(st, -6.87, 0.06)
    assert sb <= allowables.tension_full_prestress()      # 台灣零拉
    assert st >= allowables.comp_service(40)              # 頂緣壓 ≤ 0.6f'c


def test_pe_min_inverse():
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    _close(Pe_min_zero_tension(sec, ten.e, M_serv) / 1e3, 22311, 50)


def test_design_adequate():
    L = compute_losses(ten, sec, M_DC, M_DW)
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    assert L.Pe >= Pe_min_zero_tension(sec, ten.e, M_serv)     # 8組×19股（HS20）足夠


def test_shear_D1():
    """D1 腹板抗剪：以引擎算出的 Pe 串接，重現 fpc/Vp/σ1/Vcw/Vs。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    s = shear_web(Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=2069e3, x_control=1692, L=40000)
    _close(s.fpc, 4.68, 0.05)
    _close(s.Vp / 1e3, 1204, 10)         # kN（每腹板）
    _close(s.sigma1, 3.08, 0.05)
    _close(s.Vcw / 1e3, 2050, 10)        # kN
    _close(s.Vs_req / 1e3, 384, 10)      # kN
    assert not s.sigma1_ok               # 近支承主拉超限（靠箍筋）→ 設計常態


def test_shear_D16at250_passes():
    """D16@250 雙腳箍（Av/s=1.590）→ φVn=2,823 > Vu=2,329。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    s = shear_web(Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=2069e3, x_control=1692, L=40000)
    Av_s = 397.4 / 250                   # D16 雙腳 @250mm
    cap = phiVn(s.Vcw, Av_s, dv=1692)
    _close(cap / 1e3, 2702, 15)          # kN
    assert cap >= 2069e3                 # 通過
    assert Av_s >= Av_s_min_TW(40, 250)  # ≥ 最小箍筋


def test_flexure_M1():
    """M1 極限強度（8組×19股、台灣 HS20-44）：c/fps/Mn/CR。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    Mu = combinations(M_DC, M_DW, M_LL_IM)["Strength_I"]
    f = flexural_strength(ten, sec, fc=40, b_eff=8000, hf=250, dp=1880,
                          Mu_kNm=Mu, Pe=Pe, e=ten.e)
    _close(f.beta1, 0.764, 0.002)
    _close(f.c, 185, 2)
    assert f.in_flange                   # NA 在翼板內 → 矩形公式
    _close(f.fps, 1809, 5)
    _close(f.Mn, 69637, 100)             # kN·m
    _close(f.CR, 1.42, 0.02)
    assert f.phi == 1.0                  # 拉力控制
    assert f.ok                          # φMn ≥ Mu 且 ≥ 1.2Mcr 下限


def test_flexure_21strand_HL93_alt():
    """★ 引擎可處理 HL-93 重載另案（8組×21、Mu=56,856）：Mn≈76,438、CR≈1.34。
    證明同一引擎兼容兩種載重基準（台灣 HS20-44 為黃金範例，HL-93 為對照）。"""
    ten21 = Tendon(8, 21, 1109)
    Pe21 = compute_losses(ten21, sec, M_DC, M_DW).Pe
    f = flexural_strength(ten21, sec, fc=40, b_eff=8000, hf=250, dp=1880,
                          Mu_kNm=56856, Pe=Pe21, e=1109)
    _close(f.Mn, 76438, 100)
    _close(f.CR, 1.34, 0.02)
    _close(f.c, 204, 2)


def test_influence_simple_span():
    """影響線（簡支 40m）：峰值 a(L−a)/L、卡車絕對最大、車道、每車道合成。
    ★ 此為 JS 網頁計算器的共用黃金答案（peak 10m 等）。"""
    _close(il_moment_peak(40, 20), 10.0, 0.01)    # 跨中峰值
    _close(il_moment_peak(40, 10), 7.5, 0.01)
    _close(il_shear_simple(40, 40, 40), 0.0, 0.01)  # 支承 V 影響線端值
    _close(taiwan_truck_moment(40), 2860, 5)      # 台灣 HS20-44 卡車絕對最大（設計基準）
    _close(taiwan_lane_moment(40), 2680, 5)       # 車道含集中載重 80
    _close(taiwan_per_lane_moment(40), 3418, 5)   # = loads 用的每車道
    _close(hl93_per_lane_moment(40), 5673, 5)     # HL-93 對照（影響線法仍支援）


def test_taiwan_hs20_live_load():
    """台灣 HS20-44 每車道活載（卡車或車道取大，衝擊 I=15.24/(L+38.1)）。"""
    _close(taiwan_impact(40), 0.195, 0.002)
    _close(taiwan_per_lane_moment(40), 3418, 5)    # 卡車 2860 控制（< HL-93 5673）
    _close(taiwan_per_lane_shear(40), 363, 3)      # 車道 304 控制（< HL-93 588）


def test_fatigue_P1():
    """疲勞 P1：鋼腱應力幅、混凝土壓疲勞、箍筋疲勞（@250 超→@150 過）。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    fa = fatigue_check(sec, Pe, ten.e, M_perm_kNm=28800, dM_fatigue_kNm=3222, fc=40)
    _close(fa.dsig_ps, 12.6, 0.2)             # 鋼腱應力幅
    assert fa.ps_ok                            # ≤ 125
    _close(fa.sig_c_max, 6.59, 0.1)           # 混凝土壓疲勞（19股 Pe≈23,720，與原 P1 演算一致）
    assert fa.c_ok                             # ≤ 0.40f'c
    d250, ok250 = stirrup_fatigue(565, 250, 402, 1692)
    d150, ok150 = stirrup_fatigue(565, 150, 402, 1692)
    _close(d250, 208, 3)
    assert not ok250                           # @250 超限（近支承疲勞控制）
    _close(d150, 125, 3)
    assert ok150                               # @150 通過


def test_torsion_D2():
    """扭力 D2：箱梁 T_cr 極大、Tu 遠低於門檻 → 可免顯式扭設計但須閉合箍。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    tr = torsion_check(sec, Pe, fc=40, Acp=23.1e6, pc=26200, Tu_kNm=1900)
    _close(tr.fpc, 4.68, 0.05)
    _close(tr.Tcr, 42367, 100)                 # 19股 Pe → 回到 D2 演算原值 42,380
    _close(tr.threshold, 9533, 30)
    assert tr.neglect                           # Tu 1,900 << 門檻 → 可忽略
    assert tr.need_closed_stirrup               # 箱梁恆需閉合箍


def test_transverse_D3():
    """橫向 D3：懸臂/跨中/墩面 RC 板撓曲 φMn ≥ Mu。"""
    _close(slab_flexure(105.8, 1571, 200, 40, 420).phiMn, 113.0, 0.5)   # 懸臂 D20@200
    _close(slab_flexure(133.8, 2172, 200, 40, 420).phiMn, 153.2, 0.5)   # 跨中 D22@175
    _close(slab_flexure(150.3, 2534, 200, 40, 420).phiMn, 176.6, 0.5)   # 墩面 D22@150
    assert slab_flexure(150.3, 2534, 200, 40, 420).ok
    _close(As_min_slab(40, 420, 1000, 200), 571, 2)


def test_temperature_T1():
    """溫度載重定義：Zone3 正梯度 18/5 → PC 箱梁負梯度 -5.4/-1.5（×-0.30）。"""
    g = temp_gradient_AASHTO(18.0, 5.0, True)
    _close(g["neg_T1"], -5.4, 0.05)
    _close(g["neg_T2"], -1.5, 0.05)
    assert g["gamma_TG"] == 0.5


def test_bearing_E1():
    """支承 E1（HS20 純化反力）：γ_S/σ_TL限值/形狀係數/穩定/H_m/上拔 全檢核。"""
    R_LL = 290 * taiwan_per_lane_shear(40) / 588          # HS20 支承活載反力 ≈179
    b = bearing_check(1440 + R_LL, 1440, R_LL, 40, 100, 550, 450, te=10, G_kgf=8)
    _close(b.gamma_s, 0.40, 0.01)
    _close(b.shape_S, 12.4, 0.1)                          # 形狀係數
    _close(b.sigma_TL, 6.54, 0.05)                        # R_max≈1,619（純化，舊 HL-93 為 6.99）
    _close(b.sigma_TL_limit, 10.99, 0.05)                # = 112 kgf/cm²（< 1.66GS=16.2）
    _close(b.H_m, 77.7, 0.5)                              # 水平力（算例 79，捨入）
    assert b.gamma_ok and b.sigma_ok and b.stability_ok and b.H_ok and b.no_uplift


def test_anchorage_F1():
    """錨碇 F1（8組×21、Pi=32,826）：Pu/Tburst/Fspall/剝落筋。"""
    a = anchorage_check(ten.Pi/1e3, 8, 260, 2100, 4)
    _close(a.Pu, 4453, 5)
    _close(a.sum_Tburst, 3902, 10)
    _close(a.Fspall, 356, 3)
    _close(a.As_spall, 1484, 5)         # <1548 → 4-D22（19股回到原配置）
    Pult, margin, ok = spiral_local_bearing(a.Pu, 2919, 8.47, 104044, 50, 380)
    _close(Pult, 5644, 5)               # 螺旋 D16@50 局部承壓
    assert ok and margin > 1.0          # Pult > Pu（餘裕 1.27）


def test_expansion_E2():
    """伸縮縫 E2：縮短量 29.4、最大開度 49.4、Strip Seal 75。"""
    j = expansion_joint(8.8, 12.6, 8.0, 20)
    _close(j.shortening, 29.4, 0.1)
    _close(j.g_max, 49.4, 0.1)
    assert j.joint_type == "Strip Seal 75mm"


def test_temperature_SE_T1():
    """T1 自平衡應力（斷面積分）：Tu/TL、σ_SE、負梯度底板 Service +2.00 超限（控制工況）。"""
    bands = [ThermalBand(0,300,3_000_000,11.5), ThermalBand(300,400,80_000,2.5),
             ThermalBand(400,1750,1_080_000,0), ThermalBand(1750,2000,1_375_000,0)]
    fibers = [("頂板頂",0,18.0), ("底板底",2000,0.0)]
    r = self_equilibrating_stress(bands, 1.26e12, 870, 2000, fibers)
    _close(r.Tu, 6.27, 0.02)
    _close(r.TL, 39.6, 0.2)
    _close(r.sigma_pos["頂板頂"], 1.80, 0.03)
    _close(r.sigma_pos["底板底"], -5.31, 0.03)
    _close(r.sigma_neg["底板底"], 1.59, 0.03)     # 負梯度底板轉為拉
    total, ok = thermal_service_check(r.sigma_neg["底板底"], 1.2, 0.5)
    _close(total, 2.00, 0.03)
    assert not ok                                  # +2.00 > 0（自含斷面/無預力的孤島保守值）


def test_temperature_integrated_T1():
    """★ T1 接線：config A 斷面 + 引擎服務性底緣（含預力）→ 熱應力被預力吸收 → 通過。
    證明孤島 illustration 的 +2.00（自含斷面 h=2000、σ_base 無預力）為保守假象。"""
    bands = [ThermalBand(0, 250, 11000*250, 12.58), ThermalBand(250, 300, 700*50, 6.08),
             ThermalBand(300, 400, 700*100, 2.5), ThermalBand(400, 1900, 700*1500, 0),
             ThermalBand(1900, 2100, 5800*200, 0)]
    r = self_equilibrating_stress(bands, sec.I, sec.h - sec.yb, sec.h,
                                  [("底板底", 2100, 0.0)], Ec=29700)
    _close(r.TL, 14.39, 0.2)                        # 配置A 頂板薄 → TL≪自含斷面 39.6
    _close(r.sigma_neg["底板底"], 0.21, 0.03)
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    _, sb = stresses(compute_losses(ten, sec, M_DC, M_DW).Pe, sec, ten.e, M_serv)
    total, ok = thermal_service_check(r.sigma_neg["底板底"], sb, 0.5)
    _close(total, -0.80, 0.05)
    assert ok                                       # 接真實預力(底緣-0.91)後 → 全壓通過


def test_continuous_pier():
    """★ 連續梁中墩：次彎矩 M2、中墩 T 斷面 M1（NA 進腹板 c>hf → CR≈0.42 嚴重不足）。
    彙整最大缺口/真實控制工況（簡支不會出現）。"""
    M2_mid = secondary_moment(8320, primary_moment([(23700, 0.950), (12557, -0.300)]))
    M2_pier = secondary_moment(-10594, primary_moment([(23700, -0.080), (12557, 0.900)]))
    _close(M2_mid, -10428, 5)
    _close(M2_pier, -20000, 5)
    ft = flexural_strength_T(11292, 1860, 40, 1400, 200, 700, 1950, 75337)
    assert ft.flanged                          # NA 進腹板 → T 斷面
    _close(ft.c, 766, 2)
    _close(ft.fps, 1655, 5)
    _close(ft.Mn, 31888, 100)
    assert ft.CR < 0.5 and not ft.ok           # 嚴重不足（CR≈0.42）
    # 對照：簡支跨中正彎矩同公式 → 矩形(翼板內)、CR>1
    fr = flexural_strength_T(21280, 1860, 40, 8000, 250, 700, 1880, 48965)
    assert not fr.flanged and fr.CR > 1
    # B 墩服務性底緣（雙控之二）：Pe=36,257(底+頂)、e=-259(頂板PT形心上)、M_ext=-41,080
    _, sb_pier = pier_service_stress(36257e3, sec, -259, -41080)
    _close(sb_pier, -19.97, 0.05)              # 壓 19.97 > 0.45f'c=18（1.11 倍）
    assert sb_pier < -18.0                      # 超限 → B 墩底緣為控制斷面


def test_moment_envelope():
    """彎矩包絡線（簡支 40m）：峰值≈絕對最大、兩端=0、拋物線狀。"""
    env = moment_envelope_simple(40)
    peak = max(m for _, m in env)
    _close(peak, 2863, 5)                        # 峰值 ≈ 絕對最大彎矩
    _close(env[0][1], 0.0, 1)                     # 端點(a=0) M=0
    _close(env[-1][1], 0.0, 1)                    # 端點(a=L) M=0
    mid = [m for a, m in env if abs(a - 20) < 0.1][0]
    assert mid >= max(m for a, m in env if abs(a - 10) < 0.1)  # 跨中 ≥ 1/4 點


def test_deflection_C2C3():
    """撓度/預拱（8組×19、Pe≈23,720、HS20-44）：LBR 91.3%、δ_LL 11.7、預拱 ~18。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    w_LL = 56.7 * taiwan_per_lane_moment(40) / hl93_per_lane_moment(40)   # HS20 等效 ≈ 34.2
    d = deflection_analysis(40000, 29700, sec, w_DL=144, Pe=Pe, e=ten.e, w_LL=w_LL)
    _close(d.d_DL, 49.2, 0.2)
    _close(d.w_eq, 131.5, 0.5)
    _close(d.LBR * 100, 91.3, 0.3)
    _close(d.net_long_term, 12.8, 0.3)
    _close(d.d_LL, 11.7, 0.2)
    assert d.d_LL_ok                              # δ_LL ≤ L/800
    _close(d.camber, 18, 1)


def test_load_standard_drives_strand_count():
    """★ 證明載重標準決定配置：8組×19股 在 HL-93 2 車道（11,346）下底緣轉拉、
    Pe < Pe_min（故 HL-93 需增配 21 股）；台灣 HS20-44（6,837）下 19 股已全壓足夠。
    同一引擎、同一斷面，兩種載重標準給出不同的最小設計。"""
    M_HL93 = combinations(M_DC, M_DW, lane_live_load(5673, 2, 1.0))["Service_I"]   # 40,146
    L19 = compute_losses(ten, sec, M_DC, M_DW)        # ten 現為 19 股
    _, sb_hl93 = stresses(L19.Pe, sec, ten.e, M_HL93)
    assert sb_hl93 > 0                                            # HL-93 下底緣拉（19 股不足）
    assert L19.Pe < Pe_min_zero_tension(sec, ten.e, M_HL93)       # → HL-93 需 21 股
    _, sb_hs20 = stresses(L19.Pe, sec, ten.e, combinations(M_DC, M_DW, M_LL_IM)["Service_I"])
    assert sb_hs20 <= 0                                           # HS20-44 下底緣全壓（19 股足夠）


def test_tendon_profile_G1():
    """G1 鋼腱線形（8組×21股 HL-93 KB敘述軌，對齊算例_鋼腱線形設計）：
    w_eq=8Pa/L²、LBR、θ_end=4a/L、R=L²/8a、單/雙端摩擦損失。
    w_DL=8(M_DC+M_DW)/L²=144 kN/m 與參考橋自洽。"""
    ten21 = Tendon(8, 21, 1109)                          # HL-93 增配股數
    Pe21 = compute_losses(ten21, sec, M_DC, M_DW).Pe
    w_DL = 8 * (M_DC + M_DW) * 1e6 / 40000 ** 2          # N/mm = 144
    g = tendon_profile(ten21.Pi, Pe21, a=1109, L=40000, w_DL=w_DL)
    _close(g.theta_end, 0.111, 0.001)
    _close(g.R / 1000, 180.3, 0.5)                       # m
    _close(g.w_eq_transfer, 181.9, 0.5)                  # kN/m
    _close(g.w_eq_service, 141.8, 0.5)
    _close(w_DL, 144.0, 0.1)
    _close(g.LBR_transfer, 1.263, 0.01)
    _close(g.LBR_service, 0.985, 0.01)
    _close(g.fric_single_end, 0.161, 0.002)             # 單端、遠端
    _close(g.fric_dual_mid, 0.084, 0.002)               # 雙端、跨中
    assert g.R_ok                                        # R 180m >> R_min 3m


def test_stm_F2():
    """F2 STM 端橫隔版 General Zone（參考橋 8組×19股，對齊算例_STM端橫隔版設計）：
    T_burst=0.25ΣP(1−a/h)、d_burst、爆裂鋼筋 As、主壓桿 β 等級、CCT 節點承壓。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe       # ≈23,720 kN
    P_web = 14850e3                                     # 單腹板 4 組 × 3,713 kN
    r = general_zone_burst(P_web=P_web, a=1050, h=2100, e_anc=150,
                           Pe=Pe, Ac=sec.A, A_cs_strut=540000,
                           strut_force=12830e3, fc=40, fy=420)
    _close(r.sigma_pe, 4.68, 0.02)                      # Pe/Ac MPa
    _close(r.T_burst / 1e3, 1856, 5)                    # kN（單腹板）
    _close(r.d_burst, 900, 1)                           # mm
    _close(r.As_burst, 4911, 5)                         # mm²
    assert r.beta_strut_required == "confined"          # 主壓桿需 β=1.0 橫向束制
    assert r.strut_ok
    _close(f_cu(40, 0.80), 27.2, 0.1)                   # CCT 節點 f_cu
    _close(node_capacity(40, "CCT", 90000) / 1e3, 1714, 5)  # 節點 A φF_nn kN


def test_durability_N1():
    """N1 耐久性：灌漿驗收門檻（道示 17.6.6）、100年鋼筋應力限值（6.2.2/6.3.2）、
    PC疲勞 min(0.6Pu,0.75Py)、設計年限。standalone 卡，golden 取卡片 codified 限值。"""
    assert GROUT["w_c_max"] == 0.45 and GROUT["f28_min"] == 30.0
    assert GROUT["bleed_max_pct"] == 0.0 and GROUT["chloride_max_pct"] == 0.08
    assert grout_qc_check(0.40, 35, 0.0, 0.3, 0.05).all_ok            # 合格灌漿
    bad = grout_qc_check(0.50, 28, 1.0, 0.8, 0.10)
    assert (not bad.all_ok) and len(bad.failed) == 5                  # 全 5 項不合格
    assert rebar_stress_limit("常時") == 100.0
    assert rebar_stress_limit("疲勞_一般") == 180.0
    assert rebar_stress_limit("疲勞_床版翼緣") == 120.0
    _close(pc_fatigue_limit(1860, 1674), 1116.0, 0.1)                 # min(1116, 1255.5)
    assert design_life("AASHTO") == 75 and design_life("日本") == 100
    assert design_life("台灣") == (50, 100)


def test_construction_stage_H1H2():
    """H1/H2 施工階段（參考橋 8組×19股，對齊算例_40m參考橋施工階段應力歷程）：
    支架上施拉 M_sw=0 → 過平衡頂緣引張；全PT 超限、分批4組通過；脫架自重活化回壓。"""
    _close(transfer_tension_limit(32), 1.41, 0.01)
    s8 = batched_transfer(29700e3, 8, 8, sec, 1109, 32)        # S2 全 PT
    _close(s8.sigma_top, 1.86, 0.02)
    _close(s8.sigma_bot, -19.18, 0.05)
    assert (not s8.top_ok) and s8.bot_ok                       # 頂超限、底剛好過 fci32
    s4 = batched_transfer(29700e3, 4, 8, sec, 1109, 32)        # S2 分批 4 組
    _close(s4.sigma_top, 0.93, 0.02)
    _close(s4.sigma_bot, -9.59, 0.05)
    assert s4.top_ok and s4.bot_ok                             # 分批解決
    s3 = stage_stress(29700e3, sec, 1109, 24800, 32)           # S3 脫架（自重活化）
    _close(s3.sigma_top, -3.96, 0.03)
    _close(s3.sigma_bot, -9.15, 0.03)
    assert s3.sigma_top < 0 and s3.sigma_bot < 0               # 全斷面回壓


def test_cantilever_H3():
    """H3 平衡懸臂（80+80m 變深連續梁，對齊算例_懸臂工法施工階段設計）：
    變深 h(x)、逐步懸臂彎矩 Σ(G·arm)+掛籃、長期下撓 δ(1+φ)。"""
    # 變深 h(x)=2.2+2.3(x/40)²，端點驗證
    _close(variable_depth(0, 4.5, 2.2, 40), 2.20, 0.001)    # 跨中 = h_mid
    _close(variable_depth(40, 4.5, 2.2, 40), 4.50, 0.001)   # 墩 = h_pier
    _close(variable_depth(20, 4.5, 2.2, 40), 2.775, 0.001)
    # 逐步懸臂彎矩（例中力臂：自重對 x=4、掛籃對墩 CL 40.5）
    w = [643, 599, 550, 497, 439, 385, 353, 341]
    a = [x - 4 for x in [6.75, 11.25, 15.75, 20.25, 24.75, 29.25, 33.75, 38.25]]
    _close(cantilever_moment(w, a), 61661, 40)              # 自重項（算例 61,625，累積捨入）
    _close(cantilever_moment(w, a, 800, 40.5), 94061, 40)   # +掛籃（算例公布 94,025）
    # 長期下撓 δ(1+φ)
    _close(long_term_deflection(147, 2.0), 441, 1)


def test_launching_H4():
    """H4 推進 ILM（40m等跨等深，對齊算例_推進工法施工設計）：推進彎矩包絡、
    臨時置中預力(e=0)、底緣殘餘壓、束數、頂推力、滑動支承支壓。"""
    w, L, Lc = 120.0, 40.0, 14.0
    A, Zb = 4.870e6, 3.093e9
    _close(launching_cantilever_moment(w, Lc), -11760, 5)      # M⁻ 懸臂
    _close(launching_span_moment(w, L), 24000, 5)              # M⁺ 跨中
    Mpos = launching_span_moment(w, L)
    Pc = centric_prestress_required(Mpos, Zb, A, sigma_res=1.5)
    _close(Pc, 45094, 10)                                      # 算例公布 45,100
    assert n_tendons(Pc, 2510) == 18                           # ⌈45,094/2,510⌉
    _close(launching_bottom_stress(45100, A, Mpos, Zb), -1.50, 0.02)  # 底緣恰餘 1.5 壓
    _close(jacking_force(0.10, 120 * 200), 2400, 1)            # μ_s·W_total
    _close(bearing_stress(7200, 960000), 7.50, 0.02)           # 擴大後 < 8.38 ✓


def test_segmental_H5H6():
    """H5 預鑄節塊 + H6 接縫（對齊算例_預鑄節塊工法施工設計 / 算例_節塊接縫設計）：
    節塊重、拼裝期接縫壓、剪力鍵設計承載與 LS3 驗核、黏結 PT 比例。"""
    # H5（40m SBS，Ac=4.20m²）
    _close(segment_weight(4.20, 2.5, 25), 262.5, 0.5)              # 26.3t ≤ 30t
    _close(joint_min_prestress(4.20e6), 882, 1)                   # 0.21×Ac
    _close(joint_compression(4 * 480, 4.20e6), 0.457, 0.003)      # 4 束置中 PT
    # H6（80m BCM，Ac=2.85m²，道示 V_fuk·ξ 法）
    _close(joint_min_prestress(2.85e6), 598.5, 1)
    Vkey = shear_key_design_capacity(350, 0.439)                  # 台形鍵
    _close(Vkey, 153.65, 0.2)
    _close(shear_key_utilization(2850, 20, Vkey), 0.93, 0.01)     # (V_sd/20鍵)/V_key
    _close(bonded_pt_ratio(14000, 42500), 0.329, 0.003)          # ≥0.30 ✓


def test_design_inverse():
    """階段 4 反解設計庫：設計=驗算之逆，閉環對齊 40m 參考橋（純 Python 零相依）。"""
    Pe_min = Pe_min_zero_tension(sec, ten.e, combinations(M_DC, M_DW, M_LL_IM)["Service_I"])
    L = compute_losses(ten, sec, M_DC, M_DW)
    assert min_tendon_groups(Pe_min, 19, L.fpe) == 8          # 反解組數 = 實際 8 組（閉環）
    Pe21 = compute_losses(Tendon(8, 21, 1109), sec, M_DC, M_DW).Pe
    w_DL = 8 * (M_DC + M_DW) * 1e6 / 40000 ** 2
    _close(required_drape(0.985, w_DL, 40000, Pe21), 1109, 2)  # 反解垂度 = e_m（閉環）
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    Sb_min = min_section_modulus_Sb(L.Pe, sec.A, ten.e, M_serv, 0.0)
    assert sec.Sb >= Sb_min                                    # 實際斷面足夠（Sb ≥ Sb_min）
    _close(Sb_min / 1e9, 1.992, 0.01)                          # 與 Pe_min 為 σ_b=0 同式對偶


if __name__ == "__main__":
    L = compute_losses(ten, sec, M_DC, M_DW)
    c = combinations(M_DC, M_DW, M_LL_IM)
    st, sb = stresses(L.Pe, sec, ten.e, c["Service_I"])
    pem = Pe_min_zero_tension(sec, ten.e, c["Service_I"])
    print("=== 40m 參考橋黃金答案重現（台灣 HS20-44、8組×19股最小設計）===")
    print(f"  M_LL+IM(2車道) = {M_LL_IM:7.0f} kN·m   [golden 6,837]")
    print(f"  f_cgp          = {L.fcgp:7.2f} MPa    [8.60]")
    print(f"  總損失         = {L.loss_pct*100:7.1f} %      [20.1]")
    print(f"  fpe            = {L.fpe:7.0f} MPa    [1,115]")
    print(f"  Pe             = {L.Pe/1e3:7.0f} kN     [23,720]")
    print(f"  Strength I     = {c['Strength_I']:7.0f} kN·m   [48,965]")
    print(f"  Service I      = {c['Service_I']:7.0f} kN·m   [35,637]")
    print(f"  sigma_bot      = {sb:+7.2f} MPa    [-0.91]")
    print(f"  sigma_top      = {st:+7.2f} MPa    [-6.87]")
    print(f"  Pe_min(零拉)   = {pem/1e3:7.0f} kN     [22,311]")
    print(f"  設計足夠       = {'OK' if L.Pe>=pem else 'NG'}  (Pe {L.Pe/1e3:.0f} >= Pe_min {pem/1e3:.0f})")
    s = shear_web(L.Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=2329e3, x_control=1692, L=40000)
    cap = phiVn(s.Vcw, 397.4 / 250, dv=1692)
    print("\n=== D1 腹板抗剪（串接引擎算的 Pe）===")
    print(f"  fpc={s.fpc:.2f}[5.05]  Vp={s.Vp/1e3:.0f}[1,299]  sigma1={s.sigma1:.2f}[3.54] (限{s.sigma1_limit:.3f}→{'超→靠箍筋' if not s.sigma1_ok else 'OK'})")
    print(f"  Vcw={s.Vcw/1e3:.0f}[2,191]  Vs_req={s.Vs_req/1e3:.0f}[549]  phiVn(D16@250)={cap/1e3:.0f}[2,823] > Vu 2,329 -> {'OK' if cap>=2329e3 else 'NG'}")

    ten19 = Tendon(8, 19, 1109)
    L19 = compute_losses(ten19, sec, M_DC, M_DW)
    _, sb19 = stresses(L19.Pe, sec, ten19.e, c["Service_I"])
    print("\n=== ★ 改回 8組×19股 → 引擎自動證明需增配 ===")
    print(f"  損失={L19.loss_pct*100:.1f}%  Pe={L19.Pe/1e3:.0f}kN  sigma_bot={sb19:+.2f} MPa (>0=底緣拉→台灣零拉失敗)")
