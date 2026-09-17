"""連續梁中墩計算書：bridgecalc → 可列印 HTML（40+40 兩跨連續後張箱梁）。

無外部相依（純 Python + bridgecalc），可立即執行：
    python3 連續梁中墩計算書.py       # 產出 連續梁中墩計算書.html
瀏覽器開啟 → Ctrl+P → 另存 PDF。方程以 KaTeX(CDN) 渲染。

★ 2026-09-17 校正：次彎矩 M₂ 改以力法由鋼腱線形實算（原以 M_total 手填、正負號相反），
頂板腱偏心 900→646（y_t−保護層−管半徑）。結論翻轉——M₂ 在墩頂有利、跨中不利：
B 墩底緣通過，控制改為正彎矩區底緣拉應力；中墩強度仍不足（Mu 已計入 M₂）。
全部數值對齊 golden continuous_pier。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bridgecalc import (Section, flexural_strength_T, pier_service_stress,
                        parabola_seg, TendonGroup, primary_moment_at, continuous_prestress,
                        cont_tendon_segs)

# ── 40+40 兩跨連續後張箱梁（同 40m 參考斷面）──
sec = Section(A=5.065e6, I=3.287e12, yb=1329, h=2100)
E_TOP = 646.0                                   # 頂板腱墩頂偏心上限 = 771 − 75 − 50
cb = -(950 + 80) / 20 ** 2                       # 底板腱：端/墩 −80、跨中 +950（形心下為正）
bot = TendonGroup(23700, [parabola_seg(0, 40, 20, 950, cb), parabola_seg(40, 80, 60, 950, cb)])
ct = (300 + E_TOP) / 15 ** 2                     # 頂板腱：x=25~55，錨 +300、墩 −646
top = TendonGroup(12557, [parabola_seg(25, 40, 40, -E_TOP, ct), parabola_seg(40, 55, 40, -E_TOP, ct)])
G = [bot, top]
fm = continuous_prestress([40, 40], G)           # 力法
M2_pier, M2_mid, M2_15 = fm.X[1], fm.M2_at(20), fm.M2_at(15)
M1_pier, M1_mid, M1_15 = (primary_moment_at(G, x) for x in (40, 20, 15))
M2_bot = continuous_prestress([40, 40], [bot]).X[1]
M2_top = continuous_prestress([40, 40], [top]).X[1]
e15 = bot.segs[0].e(15)
e_eff = -(23700 * 80 + 12557 * E_TOP) / 36257
Mext_p, Mext_15 = -41080 + M2_pier, 14670 + 12000 + M2_15
st_p, sb_p = pier_service_stress(36257e3, sec, e_eff, Mext_p)
_, sb_p0 = pier_service_stress(36257e3, sec, e_eff, -41080)
st_15, sb_15 = pier_service_stress(23700e3, sec, e15, Mext_15)
_, sb_150 = pier_service_stress(23700e3, sec, e15, 14670 + 12000)
Mu_pier = 75337 - M2_pier                        # γ_P = 1.0（AASHTO 3.4.1）
ft = flexural_strength_T(11292, 1860, 40, 1400, 200, 700, 1975, Mu_pier)
tpa = cont_tendon_segs([40, 40], -80, 950, -600, 0.12, 0.40)   # 設計調整：底板腱墩頂抬高＋圓順
fa = continuous_prestress([40, 40], [TendonGroup(23700, tpa.segs), top])
_, sb_a15 = pier_service_stress(23700e3, sec, tpa.e_at(15), 14670 + 12000 + fa.M2_at(15))

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
<tr><td>x=15 m（0.375L）</td><td>{M1_15:,.0f}</td><td>{M2_15:+,.0f}</td><td>{M1_15+M2_15:,.0f}</td></tr>
<tr><td>跨中 x=20 m</td><td>{M1_mid:,.0f}</td><td>{M2_mid:+,.0f}</td><td>{M1_mid+M2_mid:,.0f}</td></tr>
<tr><td>B 墩 x=40 m</td><td>{M1_pier:+,.0f}</td><td>{M2_pier:+,.0f}</td><td>{M1_pier+M2_pier:+,.0f}</td></tr>
</table>
<p class="note">分項：底板腱 {M2_bot:+,.0f}（手算等效載重 w=8Pa/L²、端錨彎矩傳遞 −1/2 → 21,567 相同）、頂板腱 {M2_top:+,.0f}。<b>M₂ 全長為正彎矩：跨中不利、墩頂有利</b>——簡支梁無此項。</p>"""
sections.append(("二、連續 PT 次彎矩 M₂（力法）", sec2))

r_15 = row("正彎矩區底緣（x=15 m）", r"\sigma_b=-\frac{P_e}{A}-\frac{P_e e}{S_b}+\frac{M_{DL}+M_{LL}+M_2}{S_b}",
           f"P_e=23{{,}}700,\\ e={e15:.1f},\\ M={Mext_15:,.0f}", f"{sb_15:+.2f}{MPa}",
           "台灣完全預壓 零拉", sb_15 <= 0, "C1/連續")
sec3 = f"""{r_15}
<p class="note">🔴 <b>底緣出現 {sb_15:+.2f} MPa 拉應力</b>（AASHTO Service III 0.50√f'c=3.16 可通過）。<b>不計 M₂ 時為 {sb_150:+.2f} MPa 看似通過——漏算 M₂ 在正彎矩區偏不保守</b>。頂緣 {st_15:+.2f} MPa。</p>"""
sections.append(("三、正彎矩區服務性應力（控制項）", sec3))

r_svc = row("B 墩底緣使用階段壓應力", r"\sigma_b=-\frac{P_e}{A}-\frac{P_e e}{S_b}+\frac{M_{DL}+M_{LL}+M_2}{S_b}",
            f"P_e=36{{,}}257,\\ e={e_eff:.0f},\\ M={Mext_p:,.0f}", f"{sb_p:+.2f}{MPa}",
            "壓 ≤ 0.45f'c = 18.0 MPa", sb_p >= -18.0, "C1/連續")
sec4 = f"""{r_svc}
<p class="note">頂緣 {st_p:+.2f} MPa（無拉應力）。<b>不計 M₂ 時 σ_b={sb_p0:+.2f} MPa 超 18</b>——原算例「B 墩超限 1.11 倍」即此假象（2026-09-17 撤回）。</p>"""
sections.append(("四、B 墩服務性應力", sec4))

r_t = row("中墩負彎矩 T 斷面公稱強度", r"M_n=A_{ps}f_{ps}(d_p-\tfrac a2)+0.85f'_c(b-b_w)h_f(\tfrac a2-\tfrac{h_f}2)",
          f"c={ft.c:.0f}>h_f=200\\Rightarrow\\text{{NA進腹板}}", f"{ft.Mn:,.0f}" + BS + " kN" + BS + "cdot m",
          f"φM_n={ft.phiMn:,.0f} ≥ M_u={Mu_pier:,.0f}（CR={ft.CR:.2f}）", ft.ok, "M1/連續")
sec5 = f"""<p>M_u = 外力 75,337 − M₂ {M2_pier:,.0f}（γ_P=1.0）= {Mu_pier:,.0f} kN·m；d_p = 2,100 − (771 − {E_TOP:.0f}) = 1,975 mm。受壓區落在<b>窄底板</b>，中性軸 c={ft.c:.0f} mm &gt; h_f=200 → T 斷面。</p>
<table class="props">
<tr><td>中性軸 c</td><td>{ft.c:.0f} mm（{'T 斷面·NA 進腹板' if ft.flanged else '矩形'}）</td><td>鋼腱應力 f_ps</td><td>{ft.fps:,.0f} MPa</td></tr>
<tr><td>公稱彎矩 M_n</td><td>{ft.Mn:,.0f} kN·m</td><td>強度折減 φ</td><td>{ft.phi:.3f}（ε_t={ft.eps_t:.4f}）</td></tr>
</table>
{r_t}
<p class="note">🔴 <b>強度比 CR={ft.CR:.2f} &lt; 1 — 嚴重不足</b>：措施為增頂板 PT／配置非預力鋼筋／墩區加深或底板加寬。</p>"""
sections.append(("五、中墩負彎矩 T 斷面極限強度（控制項）", sec5))

sec6 = f"""<table class="props">
<tr><td>① 正彎矩區服務性</td><td>底緣 {sb_15:+.2f} MPa ✗（台灣零拉）</td><td>② B 墩服務性</td><td>底緣 {sb_p:+.2f} ✓</td></tr>
<tr><td>③ 中墩強度</td><td>CR={ft.CR:.2f} ✗</td><td>④ 次彎矩</td><td>M₂(B墩) {M2_pier:+,.0f}（墩頂有利、跨中不利）</td></tr>
</table>
<p class="note"><b>設計調整示範</b>：底板腱改全長分段拋物線、墩頂抬至 −600、低點 0.40L（G1 §五之補）→ M₂(B墩) {M2_pier:,.0f}→{fa.X[1]:,.0f}、x=15 m 底緣 {sb_a15:+.2f} MPa ✓。線形越接近吻合線形，M₂ 越小。中墩強度另需增頂板 PT 或非預力鋼筋。</p>"""
sections.append(("六、結論與設計措施", sec6))

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
<div class="summary">🔴 設計結論：M₂ 全長為正彎矩（B 墩 {M2_pier:+,.0f} kN·m）——<b>正彎矩區底緣 {sb_15:+.2f} MPa 拉（台灣零拉 ✗）</b>、B 墩底緣 {sb_p:+.2f} ✓、中墩強度 CR={ft.CR:.2f} ✗。
調整底板腱墩頂偏心並圓順可降 M₂；中墩須增頂板 PT／非預力鋼筋。</div>
<p class="note">本計算書全部數值由單一真理源 bridgecalc 產生並對齊 golden continuous_pier（2026-09-17 校正：M₂ 力法、頂板腱 e=646）。</p>
</body></html>"""

if __name__ == "__main__":
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "連續梁中墩計算書.html")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(OUT)
    print(f"已產生計算書：{out_path}（{len(OUT):,} bytes）")
    print(f"  M2(跨中/B墩)={M2_mid:,.0f}/{M2_pier:,.0f}  σb(x15/B墩)={sb_15:+.2f}/{sb_p:+.2f}  中墩 CR={ft.CR:.2f}")
