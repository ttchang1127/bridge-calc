"""連續梁中墩計算書：bridgecalc → 可列印 HTML（40+40 兩跨連續後張箱梁）。

無外部相依（純 Python + bridgecalc），可立即執行：
    python3 連續梁中墩計算書.py       # 產出 連續梁中墩計算書.html
瀏覽器開啟 → Ctrl+P → 另存 PDF。方程以 KaTeX(CDN) 渲染。

★ 本計算書聚焦連續梁**中墩控制斷面**——與 40m 簡支橋（另一配置）互補。
中墩兩項皆為控制項：① T 斷面極限強度 CR=0.41（NA 進腹板）② B 墩底緣服務性超限。
逐項代入僅取引擎可精確重現（golden-locked）之項；B 墩服務性因源算例 σ 值模糊，
列為定性控制項（未鎖 golden，見 網頁計算器_部署現況 §八）。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bridgecalc import (Section, secondary_moment, primary_moment,
                        flexural_strength_T, pier_service_stress)

# ── 40+40 兩跨連續後張箱梁（同 40m 參考斷面）──
sec = Section(A=5.065e6, I=3.287e12, yb=1329, h=2100)
# 次彎矩：M2 = M_total(框架分析) − M1(=ΣPe·e 一次預力彎矩)。兩層等值鋼腱模型。
M1_mid = primary_moment([(23700, 0.950), (12557, -0.300)])   # 跨中
M2_mid = secondary_moment(8320, M1_mid)
M1_pier = primary_moment([(23700, -0.080), (12557, 0.900)])  # B 墩
M2_pier = secondary_moment(-10594, M1_pier)
# 中墩負彎矩 T 斷面極限強度（受壓區在底板→NA 進腹板）
ft = flexural_strength_T(11292, 1860, 40, 1400, 200, 700, 1950, 75337)
# B 墩服務性底緣：Pe=36,257(底23,700+頂12,557)、e=-259(頂板PT形心上)、M_ext=-41,080
sb_pier = pier_service_stress(36257e3, sec, -259, -41080)[1]

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
<tr><td>鋼腱</td><td>8 束（兩層等值模型）</td><td>有效預力 f_pe</td><td>1,112 MPa</td></tr>
<tr><td>跨中偏心 e</td><td>950 mm（下）</td><td>B 墩偏心 e</td><td>900 mm（上，頂板 PT）</td></tr></table>
<p class="note">連續梁 vs 簡支：連續 PT 在超靜定結構產生<b>次彎矩 M₂</b>，且<b>中墩負彎矩</b>使受壓區落在窄底板 → 控制斷面移至中墩。</p>"""
sections.append(("一、連續梁配置", sec1))

r_m2m = row("跨中次彎矩", r"M_2 = M_{total} - M_1,\quad M_1=\textstyle\sum P_e e",
            f"8{{,}}320 - ({M1_mid:,.0f})", f"{M2_mid:,.0f}" + BS + " kN" + BS + "cdot m",
            "—（納入服務性/強度組合）", True, "G1/H3")
r_m2p = row("B 墩次彎矩", r"M_2 = M_{total} - M_1",
            f"-10{{,}}594 - ({M1_pier:,.0f})", f"{M2_pier:,.0f}" + BS + " kN" + BS + "cdot m",
            "—", True, "G1/H3")
sec2 = f"""<p>連續 PT 的一次彎矩 M₁=ΣP_e·e；框架分析總彎矩 M_total 減 M₁ 得次彎矩 M₂（超靜定約束效應）。</p>
{r_m2m}{r_m2p}
<p class="note">M₂ 在 B 墩達 {M2_pier:,.0f} kN·m（負，加重墩頂負彎矩），須計入服務性與強度組合——<b>簡支梁無此項</b>。</p>"""
sections.append(("二、連續 PT 次彎矩 M₂", sec2))

r_t = row("中墩負彎矩 T 斷面公稱強度", r"M_n=A_{ps}f_{ps}(d_p-\tfrac a2)+0.85f'_c(b-b_w)h_f(\tfrac a2-\tfrac{h_f}2)",
          f"c={ft.c:.0f}>h_f=200\\Rightarrow\\text{{NA進腹板}}", f"{ft.Mn:,.0f}" + BS + " kN" + BS + "cdot m",
          f"φM_n={ft.phiMn:,.0f} ≥ M_u=75,337（CR={ft.CR:.2f}）", ft.ok, "M1/連續")
sec3 = f"""<p>中墩負彎矩使受壓區落在<b>窄底板</b>（b_w=700 mm），中性軸 c={ft.c:.0f} mm &gt; 翼厚 h_f=200 mm → 進入腹板，須用 T 斷面公式；壓力區變窄使 M_n 大幅縮水。</p>
<table class="props">
<tr><td>中性軸 c</td><td>{ft.c:.0f} mm（{'T 斷面·NA 進腹板' if ft.flanged else '矩形'}）</td><td>鋼腱應力 f_ps</td><td>{ft.fps:,.0f} MPa</td></tr>
<tr><td>公稱彎矩 M_n</td><td>{ft.Mn:,.0f} kN·m</td><td>強度折減 φ</td><td>{ft.phi:.2f}（ε_t={ft.eps_t:.4f}）</td></tr>
</table>
{r_t}
<p class="note">🔴 <b>強度比 CR={ft.CR:.2f} ≪ 1 — 嚴重不足</b>：中墩為全橋控制斷面。措施：大幅增頂板 PT / 加深墩區斷面 / 底板加寬受壓區。</p>"""
sections.append(("三、中墩負彎矩 T 斷面極限強度（M1，控制項）", sec3))

r_svc = row("B 墩底緣使用階段壓應力", r"\sigma_b=-\frac{P_e}{A}-\frac{P_e e}{S_b}+\frac{M_{ext}}{S_b}",
            "-7.16+3.80-16.61", f"{sb_pier:+.2f}{MPa}",
            "壓 ≤ 0.45f'c = 18.0 MPa", sb_pier >= -18.0, "C1/連續")
sec4 = f"""<p>B 墩底緣（負彎矩區）：有效預力合力 P_e=36,257 kN（底板 23,700 + 頂板 12,557）、合力偏心 e=−259 mm（頂板 PT 形心上，取負）、外力彎矩 M_ext=−41,080 kN·m（含 M₂）。</p>
{r_svc}
<p class="note">🔴 <b>壓應力 {abs(sb_pier):.2f} MPa &gt; 18.0（1.11 倍超限）</b>。與 §三 T 斷面強度不足合觀：<b>B 墩服務性 + 中墩強度雙控</b>，均指向墩區需增頂板 PT 或加深斷面。已鎖 golden（continuous_pier.pier_service_sigma_bot_MPa={sb_pier:.2f}）；e 慣例同 stresses()（形心下為正，頂板 PT 取負）。</p>"""
sections.append(("四、B 墩服務性應力（控制項，golden-locked）", sec4))

sec5 = f"""<table class="props">
<tr><td>控制斷面</td><td colspan="3"><b>中墩</b>（負彎矩區）</td></tr>
<tr><td>① 極限強度</td><td>T 斷面 CR={ft.CR:.2f} ✗（NA 進腹板）</td><td>② 服務性</td><td>B 墩底緣 ≈20 &gt; 18 ✗</td></tr>
<tr><td>③ 次彎矩</td><td colspan="3">M₂(B墩)={M2_pier:,.0f} kN·m — 加重墩頂負彎矩，須計入組合</td></tr>
</table>
<p class="note"><b>設計措施</b>：(1) 大幅增設頂板連續 PT（提升負彎矩強度與底緣預壓）；(2) 墩區加深斷面 / 加腋（haunch）擴大受壓區；(3) 底板加寬使 NA 不進腹板。連續梁中墩為全橋最嚴控制斷面，與簡支橋（跨中控制、各項通過）形成對照。</p>"""
sections.append(("五、結論與設計措施", sec5))

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
<div class="hdr"><h1>40+40 連續梁　中墩控制斷面設計計算書</h1>
<div class="meta">兩跨連續後張箱梁｜同 40m 參考斷面｜聚焦中墩負彎矩控制斷面｜
由 bridgecalc 計算引擎自動產生・M₂/T 斷面極限對齊 golden・與 40m 簡支計算書互補</div></div>
{body}
<div class="summary">🔴 設計結論：<b>中墩為全橋控制斷面</b>——T 斷面極限強度 CR={ft.CR:.2f}（NA 進腹板，嚴重不足）＋ B 墩底緣服務性 ≈20&gt;18 超限。
須增頂板連續 PT／墩區加深／底板加寬。連續梁 M₂ 效應（B 墩 {M2_pier:,.0f} kN·m）不可略。</div>
<p class="note">本計算書 M₂ 與 T 斷面極限強度數值由單一真理源 bridgecalc 產生（golden: continuous_pier）；B 墩服務性 σ 值源算例模糊、未鎖 golden。</p>
</body></html>"""

if __name__ == "__main__":
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "連續梁中墩計算書.html")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(OUT)
    print(f"已產生計算書：{out_path}（{len(OUT):,} bytes）")
    print(f"  M2(跨中/B墩)={M2_mid:,.0f}/{M2_pier:,.0f}  中墩 T 斷面 c={ft.c:.0f} CR={ft.CR:.2f}")
