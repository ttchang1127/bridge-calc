"""幾何與預力配置資料類別。

單位約定（SI）：力 N、長度 mm、應力 MPa、彎矩 N·mm（介面以 kN·m 輸入，內部 ×1e6）。
對應知識庫：算例彙整_40m後張箱梁 §一（斷面）、算例_鋼腱線形設計（G1）。
"""
from dataclasses import dataclass


@dataclass
class Section:
    """箱梁斷面性質（全箱）。"""
    A: float    # 斷面積 mm²
    I: float    # 慣性矩 mm⁴
    yb: float   # 形心距底 mm（頂重斷面 yb > yt）
    h: float    # 全高 mm

    @property
    def yt(self) -> float:  # 形心距頂 mm
        return self.h - self.yb

    @property
    def St(self) -> float:  # 頂緣斷面模數 mm³
        return self.I / self.yt

    @property
    def Sb(self) -> float:  # 底緣斷面模數 mm³
        return self.I / self.yb


@dataclass
class Tendon:
    """後張鋼腱配置（全斷面）。"""
    n_tendons: int           # 鋼腱組數
    strands_per: int         # 每組股數
    e: float                 # 跨中偏心 mm（形心下為正）
    Ap_strand: float = 140   # 每股面積 mm²（15.2mm 絞線）
    fpu: float = 1860        # 極限強度 MPa
    fpj_ratio: float = 0.75  # 張拉應力比 fpj/fpu

    @property
    def Aps(self) -> float:   # 總鋼腱面積 mm²
        return self.n_tendons * self.strands_per * self.Ap_strand

    @property
    def fpj(self) -> float:   # 張拉應力 MPa
        return self.fpj_ratio * self.fpu

    @property
    def Pi(self) -> float:    # 初始張拉力 N（未扣損失）
        return self.fpj * self.Aps


# ── 箱梁翼板最小厚度（台灣 §8.9.1／§8.9.2，2026-09-20 NLM 原文）──────────────
# 「8.9.1 頂翼板 除廠製之預鑄構材及先拉法預力構件其最小厚度可減至 14cm 者外，最小頂翼板厚度
#   至少應為腹板或填角間淨距之 1/30 或 15cm。」
# 「8.9.2 底翼板 除廠製之預鑄構材及先拉法預力構件其厚度可減至 13cm 者外，最小之底板厚度
#   至少應為腹板或填角間淨距之 1/30 或 14cm。」（第八章 p.151）
# ⚠ 常見誤用：RC 箱梁底板才是「梁腹間淨距 1/16 且 ≥14cm」（第七章 §7.1.22 10.(3)b，p.113–114）；
#   **PC 箱梁頂、底板同為 1/30**。兩者混用會把 PC 底板厚度要求放大近一倍。
# ⚠ 原文「1/30 或 15cm」之「或」：本引擎取兩者大者（15cm 視為下限，因預鑄例外寫「可減至 14cm」）。
@dataclass
class SlabThicknessResult:
    clear_span: float      # 腹板（或填角）間淨距 mm
    top_req: float
    bot_req: float
    top_ok: bool
    bot_ok: bool
    rc_bot_req: float      # 對照：RC 箱梁底板 max(淨距/16, 140)
    note: str = ""


def box_slab_thickness_TW(top_t: float, bot_t: float, clear_span: float,
                          precast_or_pretensioned: bool = False) -> SlabThicknessResult:
    """台灣 §8.9.1／§8.9.2 箱梁頂、底翼板最小厚度（mm）。precast_or_pretensioned：廠製預鑄或先拉法。"""
    t_floor = 140.0 if precast_or_pretensioned else 150.0
    b_floor = 130.0 if precast_or_pretensioned else 140.0
    top_req = max(clear_span / 30.0, t_floor)
    bot_req = max(clear_span / 30.0, b_floor)
    return SlabThicknessResult(
        clear_span=clear_span, top_req=top_req, bot_req=bot_req,
        top_ok=top_t >= top_req - 1e-9, bot_ok=bot_t >= bot_req - 1e-9,
        rc_bot_req=max(clear_span / 16.0, 140.0),
        note="PC 箱梁頂底板同為淨距/30；RC 箱梁底板為淨距/16（第七章 §7.1.22）——勿混用")
