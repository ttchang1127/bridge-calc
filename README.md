# bridge-calc

橋梁影響線網頁計算器 ＋ Python 計算引擎（混合架構）。

- **網頁**（`index.html`）：影響線計算特定斷面 M/V，台灣 HS20-44＋車道載重，簡支＋連續梁，車道數選項，衝擊與多車道折減。純 JS、單檔、可離線。並連結 14 章完整設計計算書、呈現 13 項設計檢核黃金範例。
- **引擎**（`engine/bridgecalc`）：Python 權威計算核心（斷面/損失/組合/服務性 C1/剪力 D1/扭力 D2/橫向 D3/極限 M1/撓度 C2C3/疲勞 P1/支承 E1/伸縮縫 E2/錨碇 F1/溫度 T1/影響線），**回歸測試 23/23**。
- **黃金範例**：40m 後張箱梁，**台灣 HS20-44 + 8組×19股最小設計**（`make_calcbook.py` 產 14 章計算書；HL-93/21股版見知識庫算例，刻意雙軌）。
- **`golden_answers.json`**：Python 引擎與 JS 前端的**共用驗證源**——網頁開啟即顯示自驗證面板（影響線 JS 交叉驗證 + 設計檢核呈現）。

## 結構

```
index.html               # 網頁計算器（GitHub Pages 從 root serve）
golden_answers.json      # JS fetch 比對的黃金答案
engine/
  bridgecalc/            # Python 計算引擎
  tests/                 # 回歸測試
  make_golden.py         # 產 ../golden_answers.json（CI 用）
  make_calcbook.py       # 產 HTML 送審計算書
.github/workflows/ci.yml # push → 重產 golden + 回歸測試 + 新鮮度檢查
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
