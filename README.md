# bridge-calc

橋梁影響線網頁計算器 ＋ Python 計算引擎（混合架構）。

- **網頁**（`index.html`）：影響線計算特定斷面 M/V，台灣 HS20-44＋車道載重，簡支＋連續梁，車道數選項，衝擊與多車道折減。純 JS、單檔、可離線。
- **引擎**（`engine/bridgecalc`）：Python 權威計算核心（斷面/損失/組合/服務性/剪力/極限/撓度/影響線），回歸測試 13/13。
- **`golden_answers.json`**：Python 引擎與 JS 前端的**共用驗證源**——網頁開啟即顯示自驗證面板（N/N 對齊）。

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
