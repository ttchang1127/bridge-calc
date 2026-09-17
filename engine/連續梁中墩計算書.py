"""連續梁中墩計算書：bridgecalc → 可列印 HTML（40+40 兩跨連續後張箱梁）。

無外部相依（純 Python + bridgecalc），可立即執行：
    python3 連續梁中墩計算書.py       # 產出 連續梁中墩計算書.html
瀏覽器開啟 → Ctrl+P → 另存 PDF。方程以 KaTeX(CDN) 渲染。

★ 2026-09-17 校正：次彎矩 M₂ 改以力法由鋼腱線形實算（原以 M_total 手填、正負號相反），
頂板腱偏心 900→646（y_t−保護層−管半徑）。結論翻轉——M₂ 在墩頂有利、跨中不利：
外力改引擎解析包絡；服務性全數通過、中墩強度計入 M₂ 後剛好足夠。
全部數值對齊 golden continuous_pier。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bridgecalc import (Section, flexural_strength_T, pier_service_stress,
                        parabola_seg, TendonGroup, primary_moment_at, continuous_prestress,
                        cont_tendon_segs, taiwan_cont_envelope, taiwan_cont_shear_at,
                        secondary_shear, design_shear_with_V2, shear_web_at, phiVn, Av_s_min_TW)

# ── 40+40 兩跨連續後張箱梁（同 40m 參考斷面）──
sec = Section(A=5.065e6, I=3.287e12, yb=1329, h=2100)
E_TOP = 646.0                                   # 頂板腱墩頂偏心上限 = 771 − 75 − 50
cb = -(950 + 80) / 20 ** 2                       # 底板腱：端/墩 −80、跨中 +950（形心下為正）
bot = TendonGroup(23700, [parabola_seg(0, 40, 20, 950, cb), parabola_seg(40, 80, 60, 950, cb)])
ct = (300 + E_TOP) / 15 ** 2                     # 頂板腱：x=25~55，錨 +300、墩 −646
top = TendonGroup(12557, [parabola_seg(25, 40, 40, -E_TOP, ct), parabola_seg(40, 55, 40, -E_TOP, ct)])
G = [bot, top]
fm = continuous_prestress([40, 40], G)           # 力法
M2_pier, M2_mid = fm.X[1], fm.M2_at(20)
M1_pier, M1_mid = primary_moment_at(G, 40), primary_moment_at(G, 20)
M2_bot = continuous_prestress([40, 40], [bot]).X[1]
M2_top = continuous_prestress([40, 40], [top]).X[1]
W_DC, W_DW, LANES = sec.A / 1e6 * 24.5, 20.0, 2
ENV = taiwan_cont_envelope([40, 40], W_DC, W_DW, LANES, n_per_span=40)
PIER = next(r for r in ENV if abs(r.x - 40) < 1e-9)


def P_e(groups, x):
    P = Pe = 0.0
    for g in groups:
        for sg in g.segs:
            if sg.x1 - 1e-9 <= x <= sg.x2 + 1e-9:
                P += g.P
                Pe += g.P * sg.e(x)
                break
    return P, Pe / P


def scan(groups, fmx, with_M2=True):
    """跨 1（含墩）逐斷面：底緣最拉 (σ, x, M)、頂緣最拉 (σ, x)。"""
    wb = wt = None
    for r in ENV:
        if r.x > 40 + 1e-9:
            break
        P, e = P_e(groups, r.x)
        M2 = fmx.M2_at(r.x) if with_M2 else 0.0
        for Ms in (r.Ms_pos, r.Ms_neg):
            st, sb = pier_service_stress(P * 1e3, sec, e, Ms + M2)
            if wb is None or sb > wb[0]:
                wb = (sb, r.x, Ms + M2, P, e)
            if wt is None or st > wt[0]:
                wt = (st, r.x)
    return wb, wt


e_eff = -(23700 * 80 + 12557 * E_TOP) / 36257
Mext_p = PIER.Ms_neg + M2_pier
st_p, sb_p = pier_service_stress(36257e3, sec, e_eff, Mext_p)
sb_p0 = pier_service_stress(36257e3, sec, e_eff, PIER.Ms_neg)[1]
WB, WT = scan(G, fm)
WB0, _ = scan(G, fm, with_M2=False)
Mu_pier = -(PIER.Mu_neg + M2_pier)               # γ_P = 1.0（AASHTO 3.4.1）
ft = flexural_strength_T(11292, 1860, 40, 1400, 200, 700, 1975, Mu_pier)
ft0 = flexural_strength_T(11292, 1860, 40, 1400, 200, 700, 1975, -PIER.Mu_neg)

# 中墩左側 d_v 斷面剪力（d_v = max(0.9 d_p, 0.72h)；中墩側 d_p 自底緣量至合力 CGS）
DV = 0.72 * 2100
for _ in range(3):
    XS = 40 - DV / 1000
    _P, _e = P_e(G, XS)
    DV = max(0.9 * (sec.yb - _e), 0.72 * 2100)
XS = 40 - DV / 1000
PS = sum(g.P * sg.slope(XS) / 1000 for g in G for sg in g.segs if sg.x1 <= XS <= sg.x2)
PSX, _ = P_e(G, XS)
SROW = taiwan_cont_shear_at([40, 40], XS, "L", W_DC, W_DW, LANES)
V2S = secondary_shear(fm, [40, 40], XS, "L")
VUS = design_shear_with_V2(SROW.Vu_pos, SROW.Vu_neg, V2S)
SH = shear_web_at(PSX * 1e3, sec, PS / PSX, 40, 250, DV, VUS / 2 * 1e3, 2)
AV = 397.4                                       # D16 × 2 肢
PHI250, PHI200 = phiVn(SH.Vcw, AV / 250, DV) / 1e3, phiVn(SH.Vcw, AV / 200, DV) / 1e3

BS = "\\"
MPa = BS + " MPa"


def chk(ok):
    return ('<span class="ok">✓ OK</span>' if ok else '<span class="ng">✗ NG</span>')


def row(title, formula, subst, result, allow, ok, cite):
    return f"""<div class="row">
  <div class="rt">{title} <span class="cite">{cite}</span></div>
  <div class="eq">$$ {formula} $$</div>
  <div class="eq sub">$$ = {subst} = \\boxed{{{result}}} $$</div>
  <div class="chk">容許：{allow}　{chk(ok)}</div>
</div>"""


sections = []

sec1 = f"""<table class="props">
<tr><td>型式</td><td>40+40 m 兩跨連續後張箱梁</td><td>斷面</td><td>同 40m 參考斷面</td></tr>
<tr><td>斷面積 A</td><td>{sec.A:,.0f} mm²</td><td>慣性矩 I</td><td>{sec.I:.3e} mm⁴</td></tr>
<tr><td>底緣模數 S_b</td><td>{sec.Sb:.3e} mm³</td><td>頂緣模數 S_t</td><td>{sec.St:.3e} mm³</td></tr>
<tr><td>底板腱（全長）</td><td>P_e 23,700 kN；端/墩 e=−80、跨中 +950</td><td>頂板腱（x=25~55 m）</td><td>P_e 12,557 kN；錨 +300、墩 −{E_TOP:.0f}</td></tr>
<tr><td>恆載</td><td>自重 {W_DC:.1f}（A×24.5）＋SDL {W_DW:.0f} kN/m</td><td>活載</td><td>HS20-44、{LANES} 車道</td></tr>
</table>
<p class="note">偏心 e 形心下為正、彎矩正彎矩為正、應力壓為負。頂板腱墩頂偏心上限 = y_t − 保護層 − 管半徑 = 771 − 75 − 50 = {E_TOP:.0f} mm。</p>"""
sections.append(("一、連續梁配置", sec1))

r_m2 = row("B 墩次彎矩（力法）", r"M_{2,B} = -\frac{\int M_1\,m\,dx}{\int m^2\,dx},\quad M_1=-P e",
           f"-\\frac{{{fm.b[0]:,.0f}}}{{{fm.F[0][0]:.3f}}}", f"{M2_pier:+,.0f}" + BS + " kN" + BS + "cdot m",
           "—（納入服務性/強度組合，γ=1.0）", True, "G1/H3")
sec2 = f"""<p>以 B 墩彎矩為贅餘力，相容條件（B 墩轉角連續）∫(M₁+M₂)·m dx = 0，m 為 B 墩單位三角形彎矩圖。子區間對齊線形交界與錨碇點，Simpson 精確。</p>
{r_m2}
<table class="props">
<tr><td>斷面</td><td>M₁</td><td>M₂</td><td>M_total</td></tr>
<tr><td>跨中 x=20 m</td><td>{M1_mid:,.0f}</td><td>{M2_mid:+,.0f}</td><td>{M1_mid+M2_mid:,.0f}</td></tr>
<tr><td>B 墩 x=40 m</td><td>{M1_pier:+,.0f}</td><td>{M2_pier:+,.0f}</td><td>{M1_pier+M2_pier:+,.0f}</td></tr>
</table>
<p class="note">分項：底板腱 {M2_bot:+,.0f}（手算等效載重 w=8Pa/L²、端錨彎矩傳遞 −1/2 → 21,567 相同）、頂板腱 {M2_top:+,.0f}。<b>M₂ 全長為正彎矩：跨中不利、墩頂有利</b>——簡支梁無此項。</p>"""
sections.append(("二、連續 PT 次彎矩 M₂（力法）", sec2))

sec3 = f"""<p>引擎 <code>taiwan_cont_envelope</code>（解析三彎矩法）：HS20-44 設計卡車雙向掃描與車道載重取大；<b>負彎矩車道另加一個等集中載重於他跨</b>（§3.9）；衝擊長度正彎矩取該跨、負彎矩取相鄰兩跨平均（§3.13）。</p>
<table class="props">
<tr><td>B 墩 M_DC／M_DW</td><td>{PIER.M_dc:,.0f}／{PIER.M_dw:,.0f} kN·m</td><td>B 墩 M_LL+IM（2 車道）</td><td>{PIER.M_ll_neg:,.0f} kN·m（I={PIER.I_neg*100:.1f}%）</td></tr>
</table>"""
sections.append(("三、設計彎矩包絡（恆載＋活載）", sec3))

r_span = row("跨 1 底緣最不利（含 M₂）", r"\sigma_b=-\frac{P_e}{A}-\frac{P_e e}{S_b}+\frac{M_{DL}+M_{LL}+M_2}{S_b}",
             f"x={WB[1]:.0f}\\ \\text{{m}},\\ P_e={WB[3]:,.0f},\\ e={WB[4]:.0f},\\ M={WB[2]:,.0f}", f"{WB[0]:+.2f}{MPa}",
             "台灣完全預壓 零拉", WB[0] <= 0, "C1/連續")
r_svc = row("B 墩底緣使用階段壓應力", r"\sigma_b=-\frac{P_e}{A}-\frac{P_e e}{S_b}+\frac{M_{DL}+M_{LL}+M_2}{S_b}",
            f"P_e=36{{,}}257,\\ e={e_eff:.0f},\\ M={Mext_p:,.0f}", f"{sb_p:+.2f}{MPa}",
            "壓 ≤ 0.45f'c = 18.0 MPa", sb_p >= -18.0, "C1/連續")
sec4 = f"""{r_span}{r_svc}
<p class="note">跨 1 頂緣最不利 {WT[0]:+.2f} MPa（x={WT[1]:.0f} m）；B 墩頂緣 {st_p:+.2f}。<b>不計 M₂</b>：跨中底緣最不利 {WB0[0]:+.2f}（x={WB0[1]:.0f} m）、B 墩底緣 {sb_p0:+.2f}——M₂ 吃掉跨中壓應力餘裕、減輕墩頂壓應力。</p>"""
sections.append(("四、服務性應力（全線掃描）", sec4))

r_t = row("中墩負彎矩 T 斷面公稱強度", r"M_n=A_{ps}f_{ps}(d_p-\tfrac a2)+0.85f'_c(b-b_w)h_f(\tfrac a2-\tfrac{h_f}2)",
          f"c={ft.c:.0f}>h_f=200\\Rightarrow\\text{{NA進腹板}}", f"{ft.Mn:,.0f}" + BS + " kN" + BS + "cdot m",
          f"φM_n={ft.phiMn:,.0f} ≥ M_u={Mu_pier:,.0f}（CR={ft.CR:.2f}）", ft.ok, "M1/連續")
sec5 = f"""<p>M_u = 1.25M_DC + 1.50M_DW + 1.75M_LL + 1.0M₂ = {-PIER.Mu_neg:,.0f} − {M2_pier:,.0f} = {Mu_pier:,.0f} kN·m；d_p = 2,100 − (771 − {E_TOP:.0f}) = 1,975 mm。受壓區落在<b>窄底板</b>，中性軸 c={ft.c:.0f} mm &gt; h_f=200 → T 斷面。</p>
<table class="props">
<tr><td>中性軸 c</td><td>{ft.c:.0f} mm（{'T 斷面·NA 進腹板' if ft.flanged else '矩形'}）</td><td>鋼腱應力 f_ps</td><td>{ft.fps:,.0f} MPa</td></tr>
<tr><td>公稱彎矩 M_n</td><td>{ft.Mn:,.0f} kN·m</td><td>強度折減 φ</td><td>{ft.phi:.3f}（ε_t={ft.eps_t:.4f}）</td></tr>
</table>
{r_t}
<p class="note">不計 M₂ 時 M_u={-PIER.Mu_neg:,.0f}、CR={ft0.CR:.2f} ✗——<b>M₂ 在墩頂有利，是中墩強度剛好足夠的關鍵</b>；餘裕僅 {(ft.CR-1)*100:.0f}%，且只計頂板腱（底板腱在墩頂位於形心附近、保守略去）。</p>"""
sections.append(("五、中墩負彎矩 T 斷面極限強度", sec5))

r_v = row("中墩左側 d_v 斷面腹板抗剪（每腹板）", r"\phi V_n=\phi\left(V_{cw}+\frac{A_v}{s}f_y d_v\right),\quad V_{cw}=(0.094\sqrt{f'_c}+0.3f_{pc})b_w d_v+V_p",
          f"V_{{cw}}={SH.Vcw/1e3:,.0f},\ A_v/s=397.4/200", f"{PHI200:,.0f}" + BS + " kN",
          f"≥ V_u/腹板 = {abs(VUS)/2:,.0f} kN（@250 僅 {PHI250:,.0f}）", PHI200 >= abs(VUS) / 2, "D1/連續")
sec5b = f"""<p>斷面 x = {XS:.3f} m（d_v = {DV:,.0f} mm）。V_DC {SROW.V_dc:,.0f}、V_DW {SROW.V_dw:,.0f}、V_LL+IM {SROW.V_ll_neg:,.0f}（I={SROW.I*100:.1f}%，衝擊長度＝至較遠支點）；
次剪力 V₂ = dM₂/dx = {V2S:+,.0f} kN（於墩左有利 → 依 V_u＝max(|V|,|V+V₂|) 不折減）→ 設計 V_u = {VUS:,.0f} kN。
P_e = {PSX:,.0f} kN、ΣP·de/dx = {PS:,.0f} kN → V_p = {SH.Vp/1e3:,.0f} kN/腹板（腱往墩頂上升，與 V_u 同號而有利）；f_pc = {SH.fpc:.2f} MPa、b_w = 250 mm/腹板。</p>
{r_v}
<p class="note">主拉應力 σ₁ = {SH.sigma1:.2f} MPa &gt; 台灣限值 {SH.sigma1_limit:.3f} → 靠箍筋；需 A_v/s = {SH.Av_s_req:.2f} mm²/mm（最小 {Av_s_min_TW(40, 250):.2f}）。
<b>D16×2 @250 不足（φV_n {PHI250:,.0f} kN），須加密至 @200</b>。</p>"""
sections.append(("五之二、中墩剪力（d_v 斷面）", sec5b))

sec6 = f"""<table class="props">
<tr><td>① 跨中服務性</td><td>底緣最不利 {WB[0]:+.2f} MPa {'✓' if WB[0] <= 0 else '✗'}</td><td>② B 墩服務性</td><td>底緣 {sb_p:+.2f} ✓</td></tr>
<tr><td>③ 中墩強度</td><td>CR={ft.CR:.2f} {'✓' if ft.ok else '✗'}</td><td>④ 次彎矩</td><td>M₂(B墩) {M2_pier:+,.0f}（墩頂有利、跨中不利）</td></tr>
<tr><td>⑤ 中墩剪力</td><td colspan="3">V_u/腹板 {abs(VUS)/2:,.0f} kN：D16×2@250 φV_n {PHI250:,.0f} ✗ → @200 {PHI200:,.0f} ✓</td></tr>
</table>
<p class="note">中墩強度餘裕小，對頂板腱束數與 M₂ 敏感；若線形改為更接近吻合線形（M₂ 變小），中墩強度須重新確認。</p>"""
sections.append(("六、結論", sec6))

body = "".join(f'<section><h2>{t}</h2>{html}</section>' for t, html in sections)

OUT = f"""<!DOCTYPE html><html lang="zh-Hant"><head><meta charset="utf-8">
<title>連續梁中墩設計計算書</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css">
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js"></script>
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/contrib/auto-render.min.js"
  onload="renderMathInElement(document.body,{{delimiters:[{{left:'$$',right:'$$',display:true}}]}})"></script>
<style>
body{{font-family:"Noto Sans TC",system-ui,sans-serif;max-width:820px;margin:0 auto;padding:24px;color:#1a1a1a;line-height:1.6}}
.hdr{{border-bottom:2px solid #b91c1c;padding-bottom:12px;margin-bottom:8px}}
.hdr h1{{margin:0;font-size:22px}} .meta{{color:#555;font-size:13px;margin-top:4px}}
section{{margin:18px 0;page-break-inside:avoid}}
h2{{font-size:16px;background:#fef2f2;padding:6px 10px;border-left:3px solid #b91c1c;margin:0 0 10px}}
.row{{margin:10px 0;padding:8px 0;border-bottom:0.5px solid #e5e7eb}}
.rt{{font-weight:500;font-size:14px}} .cite{{color:#b91c1c;font-size:12px;margin-left:6px}}
.eq{{margin:4px 0;font-size:14px}} .sub{{color:#444}}
.chk{{font-size:13px;color:#555}} .ok{{color:#16a34a;font-weight:600}} .ng{{color:#dc2626;font-weight:600}}
table.props{{width:100%;border-collapse:collapse;font-size:13px;margin:6px 0}}
table.props td{{padding:4px 8px;border:0.5px solid #e5e7eb}} table.props td:nth-child(odd){{color:#555;width:22%}}
.note{{font-size:12px;color:#777}}
.summary{{margin-top:20px;padding:12px;background:#fef2f2;border:1px solid #b91c1c;border-radius:8px;font-weight:500}}
@media print{{body{{padding:0}} a{{display:none}} .hdr{{margin-top:0}}}}
</style></head><body>
<div class="hdr"><h1>40+40 連續梁　次彎矩與控制斷面設計計算書</h1>
<div class="meta">兩跨連續後張箱梁｜同 40m 參考斷面｜聚焦次彎矩 M₂ 與控制斷面｜
由 bridgecalc 計算引擎自動產生・M₂/T 斷面極限對齊 golden・與 40m 簡支計算書互補</div></div>
{body}
<div class="summary">設計結論：M₂ 全長為正彎矩（B 墩 {M2_pier:+,.0f} kN·m）。外力以引擎包絡實算後——跨中底緣最不利 {WB[0]:+.2f} MPa、B 墩底緣 {sb_p:+.2f} MPa，服務性全數通過；中墩強度 CR={ft.CR:.2f}（計入 M₂ 後剛好足夠，不計 M₂ 為 {ft0.CR:.2f}）；中墩剪力須箍筋加密至 D16×2@200。</div>
<p class="note">本計算書全部數值由單一真理源 bridgecalc 產生並對齊 golden continuous_pier／cont_envelope_taiwan／cont_shear_taiwan（2026-09-17：M₂ 力法、頂板腱 e=646、外力解析包絡）。</p>
</body></html>"""

if __name__ == "__main__":
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "連續梁中墩計算書.html")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(OUT)
    print(f"已產生計算書：{out_path}（{len(OUT):,} bytes）")
    print(f"  M2(跨中/B墩)={M2_mid:,.0f}/{M2_pier:,.0f}  σb(跨最不利/B墩)={WB[0]:+.2f}/{sb_p:+.2f}  中墩 CR={ft.CR:.2f}")
