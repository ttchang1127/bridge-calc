"""回歸測試：40m 參考橋（2 設計車道、8組×21股）黃金答案。

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
                        moment_envelope_simple, fatigue_check, stirrup_fatigue,
                        torsion_check)

# ── 40m 參考橋輸入 ──
sec = Section(A=5.065e6, I=3.287e12, yb=1329, h=2100)
ten = Tendon(n_tendons=8, strands_per=21, e=1109)
M_DC, M_DW = 24800, 4000
M_LL_IM = lane_live_load(per_lane_M=5673, n_lanes=2, m=1.0)   # 2 車道 → 11,346


def _close(a, b, tol):
    assert abs(a - b) <= tol, f"{a:.2f} != {b:.2f} (±{tol})"


def test_live_load():
    _close(M_LL_IM, 11346, 5)


def test_losses_nonlinear_coupling():
    L = compute_losses(ten, sec, M_DC, M_DW)
    _close(L.fcgp, 10.39, 0.05)
    _close(L.ES, 33, 1)
    _close(L.creep, 115, 2)
    _close(L.loss_pct * 100, 22.0, 0.3)
    _close(L.fpe, 1088, 3)
    _close(L.Pe / 1e3, 25590, 50)        # kN


def test_combinations():
    c = combinations(M_DC, M_DW, M_LL_IM)
    _close(c["Strength_I"], 56856, 30)
    _close(c["Service_I"], 40146, 5)
    _close(c["Service_III"], 37877, 5)


def test_service_stress():
    L = compute_losses(ten, sec, M_DC, M_DW)
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    st, sb = stresses(L.Pe, sec, ten.e, M_serv)
    _close(sb, -0.30, 0.06)              # 底緣全壓
    _close(st, -7.81, 0.06)
    assert sb <= allowables.tension_full_prestress()      # 台灣零拉
    assert st >= allowables.comp_service(40)              # 頂緣壓 ≤ 0.6f'c


def test_pe_min_inverse():
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    _close(Pe_min_zero_tension(sec, ten.e, M_serv) / 1e3, 25134, 50)


def test_design_adequate():
    L = compute_losses(ten, sec, M_DC, M_DW)
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    assert L.Pe >= Pe_min_zero_tension(sec, ten.e, M_serv)     # 8組×21股足夠


def test_shear_D1():
    """D1 腹板抗剪：以引擎算出的 Pe 串接，重現 fpc/Vp/σ1/Vcw/Vs。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    s = shear_web(Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=2329e3, x_control=1692, L=40000)
    _close(s.fpc, 5.05, 0.05)
    _close(s.Vp / 1e3, 1299, 10)         # kN（每腹板）
    _close(s.sigma1, 3.54, 0.05)
    _close(s.Vcw / 1e3, 2191, 10)        # kN
    _close(s.Vs_req / 1e3, 549, 10)      # kN
    assert not s.sigma1_ok               # 近支承主拉超限（靠箍筋）→ 設計常態


def test_shear_D16at250_passes():
    """D16@250 雙腳箍（Av/s=1.590）→ φVn=2,823 > Vu=2,329。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    s = shear_web(Pe, sec, ten.e, fc=40, bw_eff=250, dv=1692,
                  Vu=2329e3, x_control=1692, L=40000)
    Av_s = 397.4 / 250                   # D16 雙腳 @250mm
    cap = phiVn(s.Vcw, Av_s, dv=1692)
    _close(cap / 1e3, 2823, 15)          # kN
    assert cap >= 2329e3                 # 通過
    assert Av_s >= Av_s_min_TW(40, 250)  # ≥ 最小箍筋


def test_flexure_M1():
    """M1 極限強度（8組×21股、2 車道）：c/fps/Mn/CR。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    Mu = combinations(M_DC, M_DW, M_LL_IM)["Strength_I"]
    f = flexural_strength(ten, sec, fc=40, b_eff=8000, hf=250, dp=1880,
                          Mu_kNm=Mu, Pe=Pe, e=ten.e)
    _close(f.beta1, 0.764, 0.002)
    _close(f.c, 204, 2)
    assert f.in_flange                   # NA 在翼板內 → 矩形公式
    _close(f.fps, 1803, 5)
    _close(f.Mn, 76438, 100)             # kN·m
    _close(f.CR, 1.34, 0.02)
    assert f.phi == 1.0                  # 拉力控制
    assert f.ok                          # φMn ≥ Mu 且 ≥ 1.2Mcr 下限


def test_flexure_19strand_matches_legacy():
    """★ 引擎重現 M1 演算的歷史值（8組×19）：Mn≈69,650、CR≈1.23。"""
    ten19 = Tendon(8, 19, 1109)
    Pe19 = compute_losses(ten19, sec, M_DC, M_DW).Pe
    f = flexural_strength(ten19, sec, fc=40, b_eff=8000, hf=250, dp=1880,
                          Mu_kNm=56784, Pe=Pe19, e=1109)
    _close(f.Mn, 69650, 100)
    _close(f.CR, 1.23, 0.02)
    _close(f.c, 185.2, 1)


def test_influence_simple_span():
    """影響線（簡支 40m）：峰值 a(L−a)/L、卡車絕對最大、車道、每車道合成。
    ★ 此為 JS 網頁計算器的共用黃金答案（peak 10m 等）。"""
    _close(il_moment_peak(40, 20), 10.0, 0.01)    # 跨中峰值
    _close(il_moment_peak(40, 10), 7.5, 0.01)
    _close(il_shear_simple(40, 40, 40), 0.0, 0.01)  # 支承 V 影響線端值
    _close(abs_max_moment(40), 2867, 5)           # HL-93 卡車絕對最大
    _close(lane_moment_simple(40), 1860, 2)
    _close(hl93_per_lane_moment(40), 5673, 5)     # = loads.lane_live_load 的 per_lane


def test_fatigue_P1():
    """疲勞 P1：鋼腱應力幅、混凝土壓疲勞、箍筋疲勞（@250 超→@150 過）。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    fa = fatigue_check(sec, Pe, ten.e, M_perm_kNm=28800, dM_fatigue_kNm=3222, fc=40)
    _close(fa.dsig_ps, 12.6, 0.2)             # 鋼腱應力幅
    assert fa.ps_ok                            # ≤ 125
    _close(fa.sig_c_max, 6.47, 0.1)           # 混凝土壓疲勞（新 Pe；P1 演算 6.59 為舊 Pe）
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
    _close(tr.fpc, 5.05, 0.05)
    _close(tr.Tcr, 43764, 100)                 # 新 Pe（D2 演算 42,380 為舊 Pe）
    _close(tr.threshold, 9847, 30)
    assert tr.neglect                           # Tu 1,900 << 門檻 → 可忽略
    assert tr.need_closed_stirrup               # 箱梁恆需閉合箍


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
    """撓度/預拱（8組×21、Pe=25,590）：LBR 98.6%、δ_LL 19.4、預拱 ~7。"""
    Pe = compute_losses(ten, sec, M_DC, M_DW).Pe
    d = deflection_analysis(40000, 29700, sec, w_DL=144, Pe=Pe, e=ten.e, w_LL=56.7)
    _close(d.d_DL, 49.2, 0.2)
    _close(d.w_eq, 141.9, 0.5)
    _close(d.LBR * 100, 98.6, 0.3)
    _close(d.net_long_term, 2.1, 0.2)
    _close(d.d_LL, 19.4, 0.2)
    assert d.d_LL_ok                              # δ_LL ≤ L/800
    _close(d.camber, 7, 1)


def test_19strand_inadequate_shows_need_to_uprate():
    """★ 證明引擎自動捕捉非線性：改回 8組×19股 → 底緣轉拉、Pe < Pe_min。"""
    ten19 = Tendon(8, 19, 1109)
    L19 = compute_losses(ten19, sec, M_DC, M_DW)
    M_serv = combinations(M_DC, M_DW, M_LL_IM)["Service_I"]
    _, sb19 = stresses(L19.Pe, sec, ten19.e, M_serv)
    assert sb19 > 0                                            # 底緣拉（不足）
    assert L19.Pe < Pe_min_zero_tension(sec, ten19.e, M_serv)  # 需增配


if __name__ == "__main__":
    L = compute_losses(ten, sec, M_DC, M_DW)
    c = combinations(M_DC, M_DW, M_LL_IM)
    st, sb = stresses(L.Pe, sec, ten.e, c["Service_I"])
    pem = Pe_min_zero_tension(sec, ten.e, c["Service_I"])
    print("=== 40m 參考橋黃金答案重現（8組×21股、2 車道）===")
    print(f"  M_LL+IM(2車道) = {M_LL_IM:7.0f} kN·m   [golden 11,346]")
    print(f"  f_cgp          = {L.fcgp:7.2f} MPa    [10.39]")
    print(f"  總損失         = {L.loss_pct*100:7.1f} %      [22.0]")
    print(f"  fpe            = {L.fpe:7.0f} MPa    [1,088]")
    print(f"  Pe             = {L.Pe/1e3:7.0f} kN     [25,590]")
    print(f"  Strength I     = {c['Strength_I']:7.0f} kN·m   [56,856]")
    print(f"  Service I      = {c['Service_I']:7.0f} kN·m   [40,146]")
    print(f"  sigma_bot      = {sb:+7.2f} MPa    [-0.30]")
    print(f"  sigma_top      = {st:+7.2f} MPa    [-7.81]")
    print(f"  Pe_min(零拉)   = {pem/1e3:7.0f} kN     [25,134]")
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
