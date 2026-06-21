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
