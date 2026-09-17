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
                        taiwan_truck_shear, max_moment_moving,
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
                        min_tendon_groups, required_drape, min_section_modulus_Sb,
                        duct_layout, duct_spacing_required,
                        parabolic_curv_segs, friction_angle, friction_at, friction_profile,
                        tendon_forces, assign_jack,
                        parabola_seg, cont_tendon_segs, TendonGroup, primary_moment_at,
                        continuous_prestress, taiwan_cont_live_moment, cont_moment_il,
                        cont_dl_moment, taiwan_lane_reduction, pier_cap_tendon_segs,
                        top_slab_tendon_check, section_from_dims, haunch_profile, ContFlex,
                        cont_support_moments_point, taiwan_cont_envelope, cont_shear_il,
                        cont_dl_shear, taiwan_cont_live_shear, taiwan_lane_shear, shear_web_at, taiwan_cont_shear_at,
                        taiwan_rear_spacings, max_moment_moving, anchor_slip_loss,
                        pier_cap_tendon_force, gauss_integrate, friction_at, tendon_slip_loss,
                        segmented_tendon_force, segmented_friction_profile,
                        cont_shear_design_scan, groups_prestress_at, stirrup_max_spacing_TW,
                        STD_STIRRUP_SPACINGS)
from bridgecalc.influence import TW_HS20_AXLES, TW_HS20_SPACING, max_shear_moving
from bridgecalc import seismic as seis
from bridgecalc import retrofit as retro

import json
_eng = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_gp = next(p for p in (os.path.join(_eng, "golden_answers.json"),              # bridge_kb/計算引擎
                       os.path.join(os.path.dirname(_eng), "golden_answers.json"))  # bridge-calc/engine
           if os.path.exists(p))
with open(_gp, encoding="utf-8") as _f:
    golden = json.load(_f)

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
    """D1 腹板抗剪：Pe 與設計剪力都由引擎串接（V_u＝d_v 斷面實算 ÷2 腹板 ≈ 2,298）。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    Vu = taiwan_cont_shear_at([40], 1.692, "R", sec.A / 1e6 * 24.5, 20, 2).Vu_pos / 2
    _close(Vu, 2297.7, 0.5)
    s = shear_web(Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=Vu * 1e3, x_control=1692, L=40000)
    _close(s.fpc, 4.68, 0.05)
    _close(s.Vp / 1e3, 1204, 10)         # kN（每腹板）
    _close(s.sigma1, 3.57, 0.05)
    _close(s.Vcw / 1e3, 2050, 10)        # kN
    _close(s.Vs_req / 1e3, 653, 10)      # kN
    assert not s.sigma1_ok               # 近支承主拉超限（靠箍筋）→ 設計常態


def test_shear_D16at250_passes():
    """D16@250 雙腳箍（Av/s=1.590）→ φVn=2,702 > Vu=2,298（8組×19；2 車道 8組×21 為 2,823）。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    Vu = taiwan_cont_shear_at([40], 1.692, "R", sec.A / 1e6 * 24.5, 20, 2).Vu_pos / 2
    s = shear_web(Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=Vu * 1e3, x_control=1692, L=40000)
    Av_s = 397.4 / 250                   # D16 雙腳 @250mm
    cap = phiVn(s.Vcw, Av_s, dv=1692)
    _close(cap / 1e3, 2702, 15)          # kN
    assert cap >= Vu * 1e3               # 通過
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


def test_truck_both_directions():
    """軸組須雙向掃描（原向＋掉頭）。單向曾使短跨支承剪力低估至 9.8%、L/4 彎矩低估 1.8%。"""
    for L in (15, 20, 30, 38, 40, 45):
        _close(taiwan_truck_shear(L), 324 - 918 / L, 0.05)   # 後軸壓支承閉合解
    _close(taiwan_truck_shear(40), 301.1, 0.1)     # 40m：仍 < 車道 304 → 車道控制
    _close(taiwan_per_lane_shear(30), 359.1, 0.2)  # 30m：卡車控制（舊單向 326.9）
    _close(max_moment_moving(40, 10, (36.0, 144.0, 144.0), (0.0, 4.25, 8.5)), 2200.5, 0.5)  # L/4（舊單向 2160）
    _close(taiwan_truck_moment(40), 2860, 5)       # 絕對最大對稱，不受影響


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
    """★ 連續梁（算例_連續梁次彎矩，2026-09-17 校正）：M2 以力法實算、頂板腱 e=646（y_t−75−50）、
    外力以 taiwan_cont_envelope 實算。M2 全長正彎矩 → 跨中不利、墩頂有利；全線服務性通過、中墩強度 CR≈1.05。"""
    g = golden["continuous_pier"]
    cb = -(950 + 80) / 20 ** 2
    bot = TendonGroup(23700, [parabola_seg(0, 40, 20, 950, cb), parabola_seg(40, 80, 60, 950, cb)])
    ct = (300 + 646) / 15 ** 2
    top = TendonGroup(12557, [parabola_seg(25, 40, 40, -646, ct), parabola_seg(40, 55, 40, -646, ct)])
    r = continuous_prestress([40, 40], [bot, top])
    _close(r.X[1], 17354, 1)
    assert r.X[0] == 0 and r.X[2] == 0
    _close(r.M2_at(20), r.X[1] / 2, 1e-6)              # 支承間直線
    _close(primary_moment_at([bot, top], 40), 10008, 1)  # 頂板腱在形心上 → 正彎矩
    _close(primary_moment_at([bot, top], 20), -22515, 1) # 跨中無頂板腱
    # 手算等效載重（底板腱）：w=8Pa/L²=122.06 上揚 → +wL²/8；端錨 +1,896 傳遞 −1/2；減 M1
    w = 8 * 23700 * 1.030 / 40 ** 2
    hand = w * 40 ** 2 / 8 - 23700 * 0.080 / 2 - 23700 * 0.080
    _close(continuous_prestress([40, 40], [bot]).X[1], hand, 0.5)
    _close(g["M2_pier_kNm"], r.X[1], 1)
    # 外力改引擎實算（taiwan_cont_envelope）：全線服務性通過；不計 M2 在墩頂偏保守（更壓）
    assert g["pier_service_sigma_bot_MPa"] > -18.0 and g["pier_service_sigma_top_MPa"] <= 0
    assert g["pier_service_sigma_bot_noM2_MPa"] < g["pier_service_sigma_bot_MPa"]
    assert g["span_max_sigma_bot_MPa"] <= 0 and g["span_max_sigma_top_MPa"] <= 0
    assert g["span_max_sigma_bot_MPa"] > g["span_max_sigma_bot_noM2_MPa"]     # M2 吃掉跨中壓應力餘裕
    Mu = g["pier_Mu_noM2_kNm"] - r.X[1]
    _close(Mu, g["pier_Mu_kNm"], 1)
    ft = flexural_strength_T(10640, 1860, 40, 1400, 200, 700, 1975, Mu)
    assert ft.flanged                          # NA 進腹板 → T 斷面
    _close(ft.c, 716, 2)                       # A_ps＝4 束×2,660＝10,640（原誤用 P_e/f_pe 反推 11,292）
    _close(ft.Mn, 31081, 50)
    assert ft.phi == 1.0                       # ε_t 0.0053 → 拉力控制
    assert ft.ok and ft.CR < 1.1               # 計入 M2 後剛好足夠（CR≈1.03）
    assert not flexural_strength_T(10640, 1860, 40, 1400, 200, 700, 1975, g["pier_Mu_noM2_kNm"]).ok
    # 對照：簡支跨中正彎矩同公式 → 矩形(翼板內)、CR>1
    fr = flexural_strength_T(21280, 1860, 40, 8000, 250, 700, 1880, 48965)
    assert not fr.flanged and fr.CR > 1


def test_cont_envelope_taiwan():
    """連續梁解析影響線＋台灣 HS20-44：閉合解獨立驗算。"""
    g = golden["cont_envelope_taiwan"]
    L = 40.0
    # 兩等跨、單位載重於 a：M_B = −a(L²−a²)/(4L²)
    _close(cont_moment_il([L, L], 40, 20), -20 * (L * L - 400) / (4 * L * L), 1e-9)
    _close(cont_dl_moment([L, L], 130.4, 40), -130.4 * L * L / 8, 1e-6)
    _close(cont_dl_moment([L, L], 130.4, 15), 9 / 128 * 130.4 * L * L, 1e-6)
    # 墩頂負彎矩車道：均布載滿兩跨（面積 −L²/8）＋2 個 80 kN 於兩跨最負縱距 −L/(6√3)
    lv = taiwan_cont_live_moment([L, L], 40)
    _close(lv.lane_neg, 9.4 * (-L * L / 8) + 2 * 80 * (-L / (6 * 3 ** 0.5)), 0.05)
    _close(lv.I_neg, taiwan_impact(40), 1e-12)               # 相鄰兩跨平均＝40（非全長 80）
    # 單跨退化為簡支
    _close(taiwan_cont_live_moment([L], 20).lane_pos, taiwan_lane_moment(40), 1e-6)
    # 三跨：負彎矩衝擊長度＝相鄰兩跨平均 35
    _close(taiwan_cont_live_moment([30, 40, 30], 30).I_neg, taiwan_impact(35), 1e-12)
    _close(cont_dl_moment([30, 40, 30], 10, 30), cont_dl_moment([30, 40, 30], 10, 70), 1e-9)
    assert taiwan_lane_reduction(2) == 1.0 and taiwan_lane_reduction(3) == 0.9 and taiwan_lane_reduction(5) == 0.75
    _close(lv.lane_neg, g["pier_lane_neg"], 0.01)


def test_pier_cap_tendon():
    """雙系統配束：墩頂局部腱＝算例頂板腱；頂板恰容一層；與全長腱合成的 M2 可疊加。"""
    g = golden["pier_cap_tendon"]
    pc = pier_cap_tendon_segs([40, 40], 15, 300, -646)
    _close(pc.segs[0].c, (300 + 646) / 225, 1e-12)
    _close(pc.segs[0].slope(40), 0.0, 1e-12)                 # 墩心斜率 0
    _close(pc.segs[0].e(25), 300, 1e-9)
    top = TendonGroup(12557, pc.segs)
    _close(continuous_prestress([40, 40], [top]).X[1], -4213, 1)   # 同算例頂板腱
    c75 = top_slab_tendon_check(2100, 1329, 250, -646, 4, 5100, 100, 75)
    assert c75.ok and c75.e_hi == c75.e_lo == -646                # 250＝75+100+75
    assert not top_slab_tendon_check(2100, 1329, 250, -700, 4, 5100, 100, 40).in_slab
    # 線性疊加：合成 M2 ＝ 各組 M2 之和
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    full = TendonGroup(23724, tp.segs)
    _close(continuous_prestress([40, 40], [full, top]).X[1],
           continuous_prestress([40, 40], [full]).X[1] + continuous_prestress([40, 40], [top]).X[1], 1e-6)
    _close(g["dual_M2_pier_kNm"], continuous_prestress([40, 40], [full, top]).X[1], 0.1)
    assert pier_cap_tendon_segs([30, 40, 30], 20, 300, -600).lengths == [15.0, 15.0]   # ≤0.5×30


def test_variable_section_haunch():
    """變斷面：數值柔度於等斷面退化為閉合解；階梯 EI 與解析積分相符；加厚使墩頂吸引彎矩。"""
    dims = (11000, 250, 5800, 200, 350, 2, 2100)
    ref = section_from_dims(*dims)
    _close(ref.A, sec.A, 1)
    _close(ref.yb, 1329, 1)
    pr0 = haunch_profile([40, 40], *dims, 200, 8)          # 不加厚 → 走數值路徑
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    g = [TendonGroup(23724, tp.segs)]
    _close(continuous_prestress([40, 40], g, stiff=pr0).X[1], continuous_prestress([40, 40], g).X[1], 1e-6)
    f1 = ContFlex([30, 40, 30], lambda x: 1.0)
    for p in (10, 35, 77):
        for a, b in zip(f1.support_moments_point(p), cont_support_moments_point([30, 40, 30], p)):
            _close(a, b, 1e-9)
    # 階梯 EI：墩兩側 c=8 m 內 EI 加倍、兩跨均布 w=10 → X = −∫M0·m/EI ÷ ∫m²/EI（多項式解析）
    L, c, k, w = 40.0, 8.0, 2.0, 10.0
    P4 = lambda x: w / (2 * L) * (L * x ** 3 / 3 - x ** 4 / 4)
    Q3 = lambda x: x ** 3 / (3 * L * L)
    num = (P4(L - c) - P4(0)) + (P4(L) - P4(L - c)) / k
    den = (Q3(L - c) - Q3(0)) + (Q3(L) - Q3(L - c)) / k
    step = ContFlex([L, L], lambda x: k if abs(x - 40) < c else 1.0, breaks=[32, 48])
    _close(step.support_moments_dist(lambda t: w)[1], -num / den, 1e-6)
    # 加厚：形心下移、I 增、墩頂恆載彎矩增大（剛度吸引＋自重）
    pr = haunch_profile([40, 40], *dims, 400, 8)
    assert pr.dyb(40) > 0 and pr.I_rel(40) > 1 and pr.bot_t_at(36) == 300
    e_h = taiwan_cont_envelope([40, 40], 124.0925, 20, 2, n_per_span=20, stiff=pr)
    assert e_h[20].M_dc < cont_dl_moment([40, 40], 124.0925, 40)
    _close(golden["variable_section_haunch"]["haunch_env_pier_M_dc"], e_h[20].M_dc, 0.1)


def test_cont_shear_taiwan():
    """連續梁剪力：IL＝彎矩 IL 之導數、DL 閉合解、單跨退化簡支、算例墩左抗剪（D16@250 不足）。"""
    sp, h = [40, 40], 1e-5
    fd = (cont_moment_il(sp, 10 + h, 30) - cont_moment_il(sp, 10 - h, 30)) / (2 * h)
    _close(cont_shear_il(sp, 10, 30, "R"), fd, 1e-6)
    _close(cont_shear_il(sp, 40, 20, "L"), (cont_moment_il(sp, 40, 20) - cont_moment_il(sp, 40 - h, 20)) / h, 1e-5)
    _close(cont_dl_shear(sp, 10, 40, "L"), -5 * 10 * 40 / 8, 1e-9)       # 內支承左側 −5wL/8
    _close(cont_dl_shear(sp, 10, 0, "R"), 3 * 10 * 40 / 8, 1e-9)         # 端支承 3wL/8
    lvS = taiwan_cont_live_shear([40], 0, "R")
    _close(lvS.truck_pos, taiwan_truck_shear(40), 0.01)
    _close(lvS.lane_pos, taiwan_lane_shear(40), 0.01)
    # 反力＝左右剪力差：單位載重於 p=20 時 B 墩反力 = V(40R) − V(40L)
    RB = cont_shear_il(sp, 40, 20, "R") - cont_shear_il(sp, 40, 20, "L")
    _close(RB, 20 / 40 + 20 * (40 ** 2 - 20 ** 2) / (4 * 40 ** 2) * 2 / 40, 1e-9)
    g = golden["cont_shear_taiwan"]
    assert g["case_Vu_design"] == g["case_Vu_neg"]                        # V2 於墩左有利 → 不折減
    assert g["case_phiVn_D16x2_s250_kN"] < abs(g["case_Vu_design"]) / 2 < g["case_phiVn_D16x2_s200_kN"]
    assert shear_web_at(1e6, sec, -0.05, 40, 250, 1500, -1e6).Vp > 0      # 同號（墩左：腱上升、剪力負）→ 有利
    assert shear_web_at(1e6, sec, -0.05, 40, 250, 1500, 1e6).Vp < 0       # 反號 → 不利


def test_rear_axle_scan():
    """HS20-44 中後軸距 4.25～9.15：簡支恆為 4.25（與 influence.py 相同）；短跨連續梁墩頂由 9.15 控制。"""
    from bridgecalc.influence_cont import _ILCache, _truck_extremes, _GRID
    sp = taiwan_rear_spacings()
    assert sp[0] == 4.25 and sp[-1] == 9.15
    for L in (20.0, 30.0, 40.0):
        m = taiwan_cont_live_moment([L], L / 2)
        assert m.V_pos == 4.25
        _close(m.truck_pos, max_moment_moving(L, L / 2, TW_HS20_AXLES, TW_HS20_SPACING, 0.05), 1e-6)
        _close(taiwan_cont_live_shear([L], 0, "R").truck_pos, taiwan_truck_shear(L), 1e-6)
    lv = taiwan_cont_live_moment([15, 15], 15)
    c = _ILCache([15, 15])
    vals = [c.eta(15, k * _GRID) for k in range(int(round(30 / _GRID)) + 1)]
    fixed = _truck_extremes(vals, 30, (15, lambda p, sd: c.eta(15, p)), spacings=[4.25])
    assert lv.V_neg == 9.15 and lv.truck_neg < fixed[1] - 30      # 長軸距使兩後軸分落兩跨最負區


def test_simple_shear_dv():
    """簡支 d_v 斷面設計剪力：DL w(L/2−x)、車道 w(L−x)²/(2L)+P_V(L−x)/L、卡車＝max_shear_moving。"""
    g = golden["simple_shear_dv"]
    L, x = 40.0, 1.692
    w_dc = sec.A / 1e6 * 24.5
    _close(g["V_dc"], w_dc * (L / 2 - x), 1e-3)
    _close(g["lane"], 9.4 * (L - x) ** 2 / (2 * L) + 116 * (L - x) / L, 1e-3)
    _close(g["truck"], max_shear_moving(L, x), 0.5)            # 引擎 0.05 m 格點＋軸在斷面；簡支函式 0.1 m 步進
    _close(g["I"], taiwan_impact(L - x), 1e-6)                    # 衝擊長度＝至較遠支點
    ll = max(g["truck"], g["lane"]) * (1 + g["I"]) * 2
    _close(g["Vu_total"], 1.25 * g["V_dc"] + 1.5 * g["V_dw"] + 1.75 * ll, 0.01)


def test_pier_cap_tendon_losses():
    """頂板腱逐點損失：滑移面積條件、雙端對稱、錨碇處滑移最大、墩頂不受滑移；精算後 CR 仍 ≥1。"""
    g = golden["pier_cap_tendon_losses"]
    for L_m, p in ((15, 5.0), (5, 5.0), (40, 2.0)):
        r = anchor_slip_loss(0, L_m, p, 6)
        area = gauss_integrate(lambda d: anchor_slip_loss(d, L_m, p, 6).dsigma, 0, L_m, [r.L_set], 0.01)
        _close(area, 6 * 195000 / 1000, 0.5)                      # ∫Δσ = Δa·Ep
    pc = pier_cap_tendon_segs([40, 40], 15, 300, -646)
    f = lambda x: pier_cap_tendon_force(pc, x, 1395.0, 10640, g["other"])
    _close(f(30).fpe, f(50).fpe, 1e-9)                            # 雙端張拉對稱
    assert f(25).slip > 150 and f(40).slip == 0 and f(10).P == 0
    _close(f(40).friction, 1395 * friction_at(15, 30, [(0, 15, 2 * 946 / 225 / 1000), (15, 30, 2 * 946 / 225 / 1000)], 0.25, 0.003, "start"), 1e-6)
    assert 1.0 <= g["pier_CR"] < 1.06 and g["M2_pier_kNm"] > golden["continuous_pier"]["M2_pier_kNm"]


def test_cont_shear_scan():
    """全長剪力掃描：最不利＝墩側 d_v 單點檢核；分區連續涵蓋全長且對稱；間距 ≤ 允許值；錨碇前後 V_p 驟變。"""
    g, gs = golden["cont_shear_scan"], golden["cont_shear_taiwan"]
    _close(g["worst_Av_s"], gs["case_Av_s_req"], 1e-4)
    z = g["zones"]
    assert z[0][0] == 0 and z[-1][1] == 80
    assert all(abs(z[i][1] - z[i + 1][0]) < 1e-6 for i in range(len(z) - 1))
    for a, b in zip(z, reversed(z)):
        assert a[2] == b[2] and abs(a[0] - (80 - b[1])) < 0.02
    cb = -(950 + 80) / 20 ** 2
    bot = TendonGroup(23700, [parabola_seg(0, 40, 20, 950, cb), parabola_seg(40, 80, 60, 950, cb)])
    ct = (300 + 646) / 15 ** 2
    top = TendonGroup(12557, [parabola_seg(25, 40, 40, -646, ct), parabola_seg(40, 55, 40, -646, ct)])
    fm = continuous_prestress([40, 40], [bot, top])
    rows, _ = cont_shear_design_scan([40, 40], groups_prestress_at([bot, top]), fm, sec.A / 1e6 * 24.5, 20, 2,
                                     2100, 1329, 40, 250, 2, 397.4, step=2.0, sec=sec)
    for r in rows:
        assert r.s_pick is None or r.s_pick <= r.s_allow + 1e-9
        assert r.s_pick is None or r.s_pick in STD_STIRRUP_SPACINGS
        if r.halved:
            assert r.s_code_max == 300
    assert g["anchorR_Vp"] > 3 * g["anchorL_Vp"] and g["anchorL_s"] < g["anchorR_s"]
    s_max, halved, ok = stirrup_max_spacing_TW(0.34 * 40 ** 0.5 * 250 * 1512, 40, 250, 1512, 2100)
    assert halved and s_max == 300 and ok
    assert not stirrup_max_spacing_TW(0.67 * 40 ** 0.5 * 250 * 1512, 40, 250, 1512, 2100)[2]


def test_tendon_slip_full_length():
    """全長腱錨具滑移：slip_mm=0 與既有相同；雙端張拉對稱且跨中為 0（L_set<半長）；交錯張拉端部下降。"""
    g = golden["tendon_slip_G1"]
    fpj = 0.75 * 1860
    seg40 = parabolic_curv_segs(40, 1109)
    assert g["simple40_mid_unaffected"] and g["simple40_slip_mid"] == 0     # 驗證 compute_losses「跨中＝0」
    _close(tendon_slip_loss(0, 40, seg40, fpj, "both", slip_mm=6),
           tendon_slip_loss(40, 40, seg40, fpj, "both", slip_mm=6), 1e-9)   # 雙端對稱
    assert tendon_slip_loss(20, 40, seg40, fpj, "both", slip_mm=0) == 0     # 預設不計
    # 單端張拉：張拉端最大、遠端為 0（L_set < L）
    assert tendon_slip_loss(0, 40, seg40, fpj, "start", slip_mm=6) > 100
    assert tendon_slip_loss(40, 40, seg40, fpj, "start", slip_mm=6) == 0
    _close(tendon_slip_loss(1.5, 80, [(0, 80, 0.0)], fpj, "start", slip_mm=6),
           tendon_slip_loss(78.5, 80, [(0, 80, 0.0)], fpj, "end", slip_mm=6), 1e-9)
    assert g["cont80_drop_pct"] < -5 and g["cont80_Pe_pier_slip6_kN"] > g["cont80_Pe_x1_5_slip6_kN"]


def test_mid_anchor_G1():
    """中間錨碇段：分段後摩擦只在段內累積 → 80 m 超 20% 者降至 14%；段數再增受滑移限制。"""
    g = golden["mid_anchor_G1"]
    fpj = 0.75 * 1860
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    cs = [(gg.x1, gg.x2, 2 * abs(gg.c) / 1000) for gg in tp.segs]
    assert g["one_seg_max_pct"] > 20 > g["two_seg_slip6_max_pct"]          # 超門檻 → 分段後合格
    assert g["three_seg_slip6_max_pct"] > g["two_seg_slip6_max_pct"] - 0.5  # 再分段幾無改善
    # 分段後段內對稱：中間錨碇兩側等距點損失相同
    a, b = (segmented_tendon_force(x, [0, 40, 80], cs, fpj, 1.0, 0.0, 0.25, 0.003, "both", 6) for x in (30, 50))
    _close(a.friction + a.slip, b.friction + b.slip, 1e-9)
    # 單段退化：anchors=[0, L] 與不分段相同
    one = segmented_tendon_force(20, [0, 80], cs, fpj, 1.0, 0.0, 0.25, 0.003, "both", 0)
    _close(one.friction / fpj, friction_at(20, 80, cs, 0.25, 0.003, "both"), 1e-12)
    # tendon_forces(anchors=…) 與逐點一致
    ts = [{"y": 300, "x": 0, "jack": "both"}]
    tf = tendon_forces(ts, 40, 80, cs, 1329, fpj, 2660, 0.0, 0.25, 0.003, 6, anchors=[0, 40, 80])
    _close(tf.per[0]["fpe"], g["pier_fpe_two_seg"], 0.1)


def test_secondary_moments_force():
    """力法 M2：Simpson 對齊折點為精確、吻合線形 M2≡0、三跨對稱、analyzer 預設與 IL 等效載重互驗。"""
    # 吻合鋼腱：e ∝ 兩跨連續梁均布載重彎矩圖（e=120x−4x²，墩頂 −1,600）→ M2 = 0
    conc = TendonGroup(20000, [parabola_seg(0, 40, 15, 900, -4.0), parabola_seg(40, 80, 65, 900, -4.0)])
    assert abs(continuous_prestress([40, 40], [conc]).X[1]) < 1e-6
    # Simpson 精確：n_sub=1 與 20 相同
    tp = cont_tendon_segs([40, 40], 0, 1109, -600)
    grp = [TendonGroup(23724, tp.segs)]
    _close(continuous_prestress([40, 40], grp, n_sub=1).X[1], continuous_prestress([40, 40], grp, n_sub=20).X[1], 1e-6)
    g = golden["cont_tendon_force_default"]
    r = continuous_prestress([40, 40], grp)
    _close(r.X[1], g["M2_pier_kNm"], 0.1)
    _close(r.X[1], 13937 * 23724 / 23307, 5)          # analyzer IL 等效載重法（μ=K=0）
    assert [round(v, 3) for v in tp.infl] == g["infl_m"]
    # 墩頂段頂點在墩心（斜率 0）、反曲點兩側斜率相等
    for x in tp.infl:
        sl = [s.slope(x) for s in tp.segs if abs(s.x1 - x) < 1e-9 or abs(s.x2 - x) < 1e-9]
        _close(sl[0], sl[1], 1e-9)
    # 三跨對稱
    tp3 = cont_tendon_segs([30, 40, 30], 0, 900, -500)
    r3 = continuous_prestress([30, 40, 30], [TendonGroup(20000, tp3.segs)])
    _close(r3.X[1], r3.X[2], 1e-6)
    # 端錨偏心 e_end≠0：M2 端點仍為 0（力法自動計入端錨彎矩）
    tpe = cont_tendon_segs([40, 40], -200, 1109, -600)
    re = continuous_prestress([40, 40], [TendonGroup(23724, tpe.segs)])
    assert re.M2_at(0) == 0 and re.M2_at(80) == 0


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


def test_friction_profile_G1():
    """摩擦損失沿長度分佈：張拉端配置決定每個斷面的 α 與路徑長。

    交叉驗證既有 G1 值：跨中（雙端）= friction_dual_mid 0.084、
    單端遠端 = friction_single_end 0.161——同一組 μ/K/a/L 由兩條路徑算出同值。
    """
    sg = parabolic_curv_segs(40, 1109)
    _close(sg[0][2], 8 * 1109 / (40 * 40 * 1000), 1e-9)       # κ = 8a/(L²·1000)
    _close(friction_angle(sg, 20, 40), 0.1109, 1e-4)          # 跨中累積 = θ_end
    _close(friction_angle(sg, 40, 40), 0.2218, 1e-4)          # 全長 = 2θ_end
    _close(friction_at(20, 40, sg, jack="both"), 0.0840, 5e-4)
    _close(friction_at(40, 40, sg, jack="start"), 0.1609, 5e-4)

    p = {jk: friction_profile(40, sg, jack=jk) for jk in ("both", "start", "end", "alt")}
    for jk in p:                                              # 跨中四配置同值（左右對稱）
        _close(p[jk].at_mid, 0.0840, 5e-4)
    _close(p["both"].at_start, 0.0, 1e-9)                     # 雙端：兩端都是張拉端
    _close(p["both"].at_end, 0.0, 1e-9)
    _close(p["both"].avg, 0.0426, 5e-4)
    _close(p["start"].at_end, 0.1609, 5e-4)                   # 單端：遠端損失最大
    _close(p["start"].at_start, 0.0, 1e-9)
    _close(p["end"].at_start, 0.1609, 5e-4)                   # 鏡像
    _close(p["start"].avg, p["end"].avg, 1e-9)
    _close(p["alt"].avg, p["start"].avg, 1e-9)                # 交錯＝半數各走一端，平均相同
    _close(p["alt"].at_start, 0.0805, 5e-4)                   # 但端部是兩者平均，不是 0 也不是 0.161
    assert p["both"].avg < p["alt"].avg                       # 雙端張拉最省
    assert p["both"].x_max == 20.0 and p["start"].x_max == 40.0


def test_friction_into_losses():
    """fric_ratio 接進 compute_losses：不傳＝原式，傳了就用該配置的損失率。"""
    base = compute_losses(ten, sec, M_DC, M_DW)
    sg = parabolic_curv_segs(40, 1109)
    same = compute_losses(ten, sec, M_DC, M_DW,
                          fric_ratio=friction_at(20, 40, sg, jack="both"))
    _close(same.Pe, base.Pe, 3e3)                             # 跨中雙端 ≡ 預設 α/x_ctrl 式
    far = compute_losses(ten, sec, M_DC, M_DW,
                         fric_ratio=friction_at(40, 40, sg, jack="start"))
    assert far.Pe < base.Pe                                   # 單端遠端損失大 → Pe 低
    _close(far.friction - base.friction,
           ten.fpj * (0.1609 - 0.0840), 2.0)                  # 差額＝損失率差×fpj


def test_tendon_forces_G1():
    """逐腱合成：總 Pe 與平均損失率等價，差別在合力位置。

    交錯端張拉時各腱損失不同：若交錯落在不同層（單列多層 → 序號＝層號），
    低層自起點、損失小、力大 → 合力下移；同層左右配對（多列）則垂直相消、
    只剩橫向偏心。全部同端張拉時各腱損失相同 → 兩者皆為 0。
    """
    sg = parabolic_curv_segs(40, 1109)

    def tf(rows, mode):
        jk = assign_jack(len(rows), mode)
        ts = [{"no": "T%d" % (i + 1), "y": rows[i][0], "x": rows[i][1], "jack": jk[i]}
              for i in range(len(rows))]
        return tendon_forces(ts, 10, 40, sg, sec.yb, 1395, 140 * 19, 180.0)

    A = [(-80, 0), (120, 0), (320, 0), (520, 0)]           # 單列四層
    B = [(150, -85), (150, 85), (290, -85), (290, 85)]     # 兩列兩層
    a, b, st = tf(A, "alt"), tf(B, "alt"), tf(A, "start")

    _close(a.Pe_total, a.Pe_avg_ratio, 1.0)                # 總力與平均損失率等價
    _close(b.Pe_total, b.Pe_avg_ratio, 1.0)
    _close(a.de, 5.10, 0.05)                               # 單列交錯 → 合力下移
    _close(a.x_eff, 0.0, 1e-9)
    _close(b.de, 0.0, 1e-9)                                # 同層配對 → 垂直相消
    _close(b.x_eff, -4.34, 0.05)                           # 但橫向偏心
    _close(st.de, 0.0, 1e-9)                               # 同端張拉 → 兩者皆 0
    _close(st.x_eff, 0.0, 1e-9)
    assert st.Pe_total > a.Pe_total                        # 全自起點在 L/4 損失較小

    assert assign_jack(4, "alt") == ["start", "end", "start", "end"]
    assert assign_jack(3, "both") == ["both"] * 3
    assert len(a.per) == 4 and a.per[0]["jack"] == "start"
    _close(a.per[0]["ratio"], friction_at(10, 40, sg, jack="start"), 1e-12)


def test_duct_layout_G1():
    """G1 管道實配：參考橋 8 組 / 2 腹板 / φ100 / 腹板 350。

    淨距規則直接決定排得下幾列，兩案並列：
    ・rule="tw"（台灣 §8.25.2 明文 max(40,1.5d_agg)=40，保護層 §8.25.1 的 40）
      → 可用寬 270、兩列兩層，e=1109 的最底層管外緣還有 100 mm ✓
    ・rule="od"（淨距再取 ≥孔道外徑=100，保護層取 PTBG 實務 75）
      → 可用寬 200、單列四層，同一配置的最底層管外緣落在梁底以下 130 mm ✗
    """
    d = duct_layout(8, 2, sec.yb - 1109, duct_od=100, web_t=350, cover=40,
                    s_v=40, d_agg=25, y_b=sec.yb, h=2100, rule="tw")
    assert duct_spacing_required(100, 25) == 40             # 預設 tw：骨材與 40 mm 控制
    assert duct_spacing_required(100, 25, "od") == 100      # od：孔道外徑控制
    assert (d.n_col, d.n_row) == (2, 2)
    _close(d.web_avail, 270, 0.1)
    _close(d.s_h, 70, 0.1)
    _close(d.pitch_v, 140, 0.1)
    _close(d.y_bot, 150, 0.1)
    _close(d.cover_bot, 100, 0.1)
    assert d.cover_ok and d.top_ok and d.fits
    _close(d.e_max, 1169, 0.5)

    ds = duct_layout(8, 2, sec.yb - 1109, duct_od=100, web_t=350, cover=75,
                     s_v=100, d_agg=25, y_b=sec.yb, h=2100, rule="od")
    assert (ds.n_col, ds.n_row) == (1, 4)
    _close(ds.cover_bot, -130, 0.1)
    assert not ds.cover_ok and not ds.fits
    _close(ds.e_max, 904, 0.5)
    _close(ds.e_max_point, 1204, 0.5)                       # 單點簡化上限（與規則無關）
    assert ds.e_max < ds.e_max_point


def test_duct_layout_pier_section():
    """墩頂斷面：腱在形心上方，控制的是頂緣保護層而非底緣。

    e_pier = −900（形心上方 900）使頂層管外緣穿出梁頂；可行下限由
    y_cgs ≤ h − cover − od/2 − (y_top − y_cgs) 反推。
    """
    d = duct_layout(8, 2, sec.yb + 900, duct_od=100, web_t=350, cover=40,
                    s_v=40, d_agg=25, y_b=sec.yb, h=2100, rule="tw")
    assert not d.top_ok and not d.fits                      # 穿出頂緣
    assert d.cover_ok                                       # 底緣反而很寬裕
    lim = sec.yb - 2100 + 40 + 50 + (d.y_top - d.y_cgs)
    _close(lim, -611, 0.5)
    ok = duct_layout(8, 2, sec.yb - lim, duct_od=100, web_t=350, cover=40,
                     s_v=40, d_agg=25, y_b=sec.yb, h=2100, rule="tw")
    _close(2100 - (ok.y_top + 50), 40, 0.5)                 # 恰落在頂緣保護層上
    assert ok.fits


def test_duct_layout_centroid_preserved():
    """實配形心必須恰等於分析用的 CGS——3D 畫多腱不得改變引擎在用的偏心 e。

    含餘數層（5 腱 / 2 列 → 2+2+1）；不對稱填充時平移量不是層距的整數倍。
    採 rule="od" 使 φ90 的淨距需求為 90，才排得出兩列的餘數情形。
    """
    for n, nw, wt, sv in ((8, 2, 350, 100), (5, 1, 500, 100), (12, 3, 600, 60), (7, 2, 450, 50)):
        d = duct_layout(n, nw, 260, duct_od=90, web_t=wt, cover=50, s_v=sv, d_agg=20, rule="od")
        ys = [y for y, c in d.rows for _ in range(c)]
        assert len(ys) == d.n_per_web
        _close(sum(ys) / len(ys), 260, 1e-6)              # 形心守恆
        assert d.rows[0][0] == min(ys) and d.rows[-1][0] == max(ys)
    d5 = duct_layout(5, 1, 260, duct_od=90, web_t=500, cover=50, s_v=100, d_agg=20, rule="od")
    assert [c for _, c in d5.rows] == [2, 2, 1]           # 由下往上逐層填滿


def test_duct_layout_checks():
    """三項構造檢核各自可獨立觸發，且 e_max 反推自洽（採 rule="od" 保守側，需求 100）。"""
    base = dict(duct_od=100, web_t=500, cover=75, d_agg=25, y_b=1329, h=2100, rule="od")
    wide = duct_layout(8, 2, 1329 - 1109, s_v=100, **base)
    assert wide.n_col == 2 and wide.n_row == 2            # 腹板加厚 → 兩列兩層
    _close(wide.s_h, 150, 0.1)

    loose = duct_layout(8, 2, 1329 - 1109, s_v=40, **base)
    assert not loose.s_v_ok                               # 40 < 需求 100
    assert loose.cover_ok                                 # 但層距變小 → 保護層反而過

    thin = duct_layout(8, 2, 1329 - 1109, s_v=100, duct_od=100, web_t=220,
                       cover=75, d_agg=25, y_b=1329, h=2100, rule="od")
    assert not thin.s_h_ok                                # 可用寬 70 < 管徑 100

    # e = e_max 時，最底層管外緣恰落在保護層上
    at_max = duct_layout(8, 2, 1329 - wide.e_max, s_v=100, **base)
    _close(at_max.cover_bot, 75, 0.5)

    tall = duct_layout(8, 1, 1329, s_v=100, duct_od=100, web_t=250,
                       cover=75, d_agg=25, y_b=1329, h=2100, rule="od")
    assert (tall.n_col, tall.n_row) == (1, 8)             # 薄腹板單列 → 8 腱疊 8 層
    assert not tall.top_ok                                # 疊高 1,400 → 頂層穿出斷面頂緣


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
    assert (not s8.top_ok) and (not s8.bot_ok)                 # 頂拉超限、底壓亦超限 0.55f'ci=17.6（裁示0.55）
    s4 = batched_transfer(29700e3, 4, 8, sec, 1109, 32)        # S2 分批 4 組
    _close(s4.sigma_top, 0.93, 0.02)
    _close(s4.sigma_bot, -9.59, 0.05)
    assert s4.top_ok and s4.bot_ok                             # 分批解決（底 −9.59 < 17.6 ✓）
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
    _close(cantilever_moment(w, a, 800, 36.5), 90861, 40)   # +掛籃（力臂統一對 x=4；算例 90,825）
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


def test_seismic_S1_falloff():
    """S1 落橋防止：min N_L(式8-10)、u_G、L_N 需求。對齊算例_落橋防止系統設計。"""
    _close(seis.min_falloff_length(40, 10, 0), 70.0, 0.01)      # (50+10+10)·1 = 70cm
    _close(seis.min_falloff_length(40, 10, 30), 77.88, 0.01)    # 斜角修正 (1+30²/8000)
    _close(seis.ground_relative_displacement("第二類", 5000, 1.2), 22.5, 0.01)  # ε_G·Le·ratio
    _close(seis.required_falloff_length(70, 15, 22.5), 70.0, 0.01)  # max(70,37.5)=min N_L 控制
    _close(seis.restrainer_yield_strength(1800), 2700.0, 0.1)   # F_y=1.5·R_d


def test_seismic_S2_isolation():
    """S2 隔震等效線性化迭代收斂 + 表3-1 阻尼修正內插。對齊算例_隔震與消能設計。"""
    _close(seis.damping_correction_B1(0.15), 1.375, 1e-4)       # 表3-1 線性內插
    _close(seis.damping_correction_BS(0.15), 1.465, 1e-4)
    r = seis.isolation_design(8000, 400, 6000, 0.60)            # 剛性墩、一般工址
    _close(r.D_d * 1000, 221.2, 0.5)                            # 收斂設計位移
    _close(r.T_e, 2.030, 0.005)                                 # 有效週期
    _close(r.xi_e * 100, 14.75, 0.1)                            # 系統等效阻尼
    _close(r.V_b_secant, 1727, 2)                               # 式C7-3 = 式C7-4（剛性墩）
    _close(r.V_b_secant, r.V_b_bilinear, 1)
    assert r.iterations <= 10                                   # 快速收斂


def test_seismic_S3_ductility():
    """S3 橋墩韌性容量設計：M_p、V_u、圍束 ρ_s(式5-5/5-6)、間距/塑鉸長度。對齊算例_橋墩韌性耐震設計。"""
    _close(seis.overstrength_moment(3000), 3900, 0.1)           # M_p=1.3Mn
    _close(seis.capacity_shear(3900, 8), 487.5, 0.1)            # V_u=M_p/H
    _close(seis.rho_s_circular(280, 2800, 17671, 15394, 800000) * 100, 0.843, 0.005)  # 式5-6控制
    # 式5-6（軸力式）> 式5-5（幾何式）
    assert 0.12 * (280/2800) * (0.5 + 1.25*800000/(280*17671)) > 0.45 * (280/2800) * (17671/15394 - 1)
    _close(seis.confinement_spacing_limit(150, 3.6), 15.0, 0.01)   # min(15,37.5,21.6)
    _close(seis.plastic_hinge_length(150, 800), 150.0, 0.01)      # max(150,133,45)


def test_seismic_S5_liquefaction():
    """S5 液狀化土壤參數折減 D_E（表8-1 三級距×深度×R_s）。對齊公式卡_液狀化與基礎耐震。"""
    _close(seis.liquefaction_reduction_DE(0.3, 5, 0.2), 0.0, 1e-9)     # 第一級淺層鬆砂→參數設零
    _close(seis.liquefaction_reduction_DE(0.3, 5, 0.4), 1/6, 1e-6)     # 密砂折減較輕
    _close(seis.liquefaction_reduction_DE(0.3, 15, 0.2), 1/3, 1e-6)    # 深層(>10m)
    _close(seis.liquefaction_reduction_DE(0.5, 5, 0.2), 1/3, 1e-6)     # 第二級
    _close(seis.liquefaction_reduction_DE(0.8, 5, 0.4), 1.0, 1e-9)     # 第三級密砂不折減
    _close(seis.liquefaction_reduction_DE(1.2, 5, 0.2), 1.0, 1e-9)     # F_L≥1 不液化


def test_retrofit_shared_cracked():
    """R2/R4 同原梁 → 共用彈性開裂換算斷面 x₁（非極限應力塊）。JTG 二次受力核心。"""
    a = 6.667
    x1 = retro.cracked_na_depth(400, 750, 1964, a)
    _close(x1, 191.3, 0.2)                                     # √(A1²+B1)−A1
    Icr = retro.cracked_inertia(400, x1, 1964, 750, a)
    _close(Icr / 1e9, 5.021, 0.02)
    _close(retro.initial_concrete_strain(200, x1, 3.0e4, Icr), 2.54e-4, 2e-6)


def test_retrofit_R1_cfrp():
    """R1 碳纖維CFRP 抗彎（案②，式6-42）。對齊算例_碳纖維CFRP抗彎補強設計。"""
    epsf = retro.cfrp_allowable_strain(2, 2.4e5, 0.167, 0.0155)
    _close(epsf, 0.007, 1e-6)                                  # 0.007 絕對上限控制
    _close(retro.cfrp_km1(2, 2.4e5, 0.167), 0.8127, 1e-3)
    r1 = retro.cfrp_moment_capacity(400, 800, 750, 13.8, 0.0033, 1964, 330,
                                    100.2, 2.4e5, epsf, 0.0003)
    _close(r1.xi_fb, 0.249, 0.001)
    _close(r1.x, 147.9, 0.3)
    assert r1.case2                                            # x ≤ ξ_fb·h → 案②
    _close(r1.Mu_kNm, 539.4, 0.5)                              # +20%
    # n_f=1 不足（<520 需求）→ 須 2 層
    epsf1 = retro.cfrp_allowable_strain(1, 2.4e5, 0.167, 0.0155)
    r1a = retro.cfrp_moment_capacity(400, 800, 750, 13.8, 0.0033, 1964, 330,
                                     50.1, 2.4e5, epsf1, 0.0003)
    assert r1a.Mu_kNm < 520


def test_retrofit_R2_plate():
    """R2 外貼鋼板 抗彎（式6-26，鋼板降伏）。對齊算例_外貼鋼板抗彎補強設計。"""
    a = 6.667
    x1 = retro.cracked_na_depth(400, 750, 1964, a)
    Icr = retro.cracked_inertia(400, x1, 1964, 750, a)
    ec1 = retro.initial_concrete_strain(200, x1, 3.0e4, Icr)
    r2 = retro.plate_moment_capacity(400, 800, 750, 13.8, 0.0033, 1964, 330,
                                     800, 305, 2.06e5, x1, ec1)
    assert r2.plate_yields                                     # ε_sp 需求 ≫ f_sp
    _close(r2.sigma_sp, 305, 0.5)
    _close(r2.x, 161.6, 0.3)
    _close(r2.Mu_kNm, 609, 1)                                  # +36%（近 40% 上限）
    _close(retro.plate_dev_length(305, 800, 2.5, 200), 788, 1)  # l_p 式6-37


def test_retrofit_R4_enlargement():
    """R4 增大截面 抗彎（式6-2/6-3）。對齊算例_增大斷面補強設計。"""
    a = 6.667
    x1 = retro.cracked_na_depth(400, 750, 1964, a)
    Icr = retro.cracked_inertia(400, x1, 1964, 750, a)
    ec1 = retro.initial_concrete_strain(200, x1, 3.0e4, Icr)
    h0 = (1964 * 750 + 982 * 860) / 2946
    _close(h0, 786.7, 0.2)
    r4 = retro.enlargement_moment_capacity(400, 860, h0, 13.8, 0.0033, 1964, 330,
                                           982, 330, 2.0e5, x1, ec1)
    assert r4.added_bar_yields
    _close(r4.x, 176.1, 0.3)
    _close(r4.Mu_kNm, 679, 1)                                  # +52%（超鋼板 40% 上限）


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
    Vu_dv = taiwan_cont_shear_at([40], 1.692, "R", sec.A / 1e6 * 24.5, 20, 2).Vu_pos / 2
    s = shear_web(L.Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=Vu_dv * 1e3, x_control=1692, L=40000)
    cap = phiVn(s.Vcw, 397.4 / 250, dv=1692)
    print("\n=== D1 腹板抗剪（串接引擎算的 Pe 與 V_u）===")
    print(f"  Vu(d_v)/腹板={Vu_dv:.0f}[2,298]")
    print(f"  fpc={s.fpc:.2f}[5.05]  Vp={s.Vp/1e3:.0f}[1,299]  sigma1={s.sigma1:.2f}[3.46] (限{s.sigma1_limit:.3f}→{'超→靠箍筋' if not s.sigma1_ok else 'OK'})")
    print(f"  Vcw={s.Vcw/1e3:.0f}[2,192]  Vs_req={s.Vs_req/1e3:.0f}[512]  phiVn(D16@250)={cap/1e3:.0f}[2,823] > Vu -> {'OK' if cap>=Vu_dv*1e3 else 'NG'}")

    ten19 = Tendon(8, 19, 1109)
    L19 = compute_losses(ten19, sec, M_DC, M_DW)
    _, sb19 = stresses(L19.Pe, sec, ten19.e, c["Service_I"])
    print("\n=== ★ 改回 8組×19股 → 引擎自動證明需增配 ===")
    print(f"  損失={L19.loss_pct*100:.1f}%  Pe={L19.Pe/1e3:.0f}kN  sigma_bot={sb19:+.2f} MPa (>0=底緣拉→台灣零拉失敗)")
