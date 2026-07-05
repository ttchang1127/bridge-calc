# bridge-calc

橋梁影響線網頁計算器 ＋ Python 計算引擎（混合架構）。

- **網頁**（`index.html`）：影響線計算器（台灣 HS20-44＋車道，簡支＋連續梁）＋**互動箱梁分析器**（輸入尺寸/鋼腱/載重 → 16 項檢核即時重算＋參數化 SVG）。純 JS、單檔、可離線。連結 16 章完整設計計算書與連續梁中墩計算書。
- **引擎**（`engine/bridgecalc`）：Python 權威計算核心 **26 模組・零相依**——箱梁核心（斷面/損失/組合/服務性 C1/剪力 D1/扭力 D2/橫向 D3/極限 M1/撓度 C2C3/疲勞 P1/支承 E1/伸縮縫 E2/錨碇 F1/溫度 T1/影響線）＋連續梁中墩＋施工 H1-H6＋反解設計庫＋**耐震 seismic（S1/S2/S3/S5，台灣單軌）**＋**補強 retrofit（R1/R2/R4，中國 JTG 單軌）**＋容許值 SSOT。**回歸 41/41、JS 38/38**。
- **黃金範例**：40m 後張箱梁，**台灣 HS20-44 + 8組×19股最小設計**（`make_calcbook.py` 產 16 章計算書＋`連續梁中墩計算書.py`；HL-93/21股版見知識庫算例，刻意雙軌）。
- **`golden_answers.json`**：Python 引擎與 JS 前端的**共用驗證源**——網頁開啟即顯示自驗證面板（影響線 JS 交叉驗證 + 設計檢核呈現）。

## 結構

```
index.html               # 網頁計算器（GitHub Pages 從 root serve）
golden_answers.json      # JS fetch 比對的黃金答案
engine/
  bridgecalc/            # Python 計算引擎（26 模組，含 seismic / retrofit）
  tests/                 # 回歸測試（41/41）
  make_golden.py         # 產 ../golden_answers.json（CI 用，勿手改 golden）
  make_calcbook.py       # 產 16 章 HTML 送審計算書
  連續梁中墩計算書.py     # 產連續梁中墩計算書 HTML
box-girder.html          # 互動箱梁分析器（Alpine，16 檢核即時重算）
box-girder-engine.js     # JS 引擎（對 golden 38/38）
.github/workflows/ci.yml # push → 重產 golden + pytest + JS + 三新鮮度檢查
```

## 本機使用

```bash
cd engine
python3 tests/test_reference_bridge.py   # 回歸測試（無需安裝）
python3 make_golden.py                    # 重產 golden_answers.json
python3 make_calcbook.py                  # 產計算書 HTML
```
網頁：用瀏覽器開 `index.html`（自驗證面板需經 http 才能 fetch；本機可用 `python3 -m http.server` 後開 localhost）。

## 部署（GitHub Pages）

Settings → Pages → Deploy from branch `main` `/ (root)`。約 1 分鐘後即可由 `https://<帳號>.github.io/bridge-calc/` 存取。

> 詳見知識庫 `部署SOP_GitHub.md`、`網頁計算器SOP.md`。核心零安裝（Python 內建 + vanilla JS 免 build）。
