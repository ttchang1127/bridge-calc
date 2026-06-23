"""溫度梯度（T1）：載重定義 + 自平衡應力斷面積分。

對應公式卡 T1、算例_溫度梯度應力設計。AASHTO Art. 3.12.3。

自平衡應力 σ_SE：非線性溫度剖面強迫斷面維持平面 → 斷面自我約束。
作法：分層積分求等效均勻溫度 Tu 與等效線性溫差 TL，
T_SE(y)=T(y)−[Tu+TL·(ȳ_t−y)/h]，σ_SE(y)=−E_c·α·T_SE(y)（+ 為拉）。
靜定梁中 Tu/TL 不生應力（自由伸縮+翹曲），故 σ_TG = σ_SE。

⚠️ T1 算例用自含斷面（h=2,000、f'c=35），非 40m 參考橋基準配置 A。
"""
from dataclasses import dataclass


def temp_gradient_AASHTO(zone_T1: float, zone_T2: float, pc_box: bool = True):
    """AASHTO 溫度梯度（§3.12.3）。

    正梯度由氣候區給 T1（頂面）、T2（300mm 深）。
    負梯度 = 正梯度 × 係數（PC 箱梁 −0.30；T 梁 −0.50）。
    回傳 dict。γ_TG（服務性溫度載重因子）= 0.5（FHWA）。
    """
    neg = -0.30 if pc_box else -0.50
    return {
        "pos_T1": zone_T1, "pos_T2": zone_T2,
        "neg_factor": neg,
        "neg_T1": round(neg * zone_T1, 2), "neg_T2": round(neg * zone_T2, 2),
        "gamma_TG": 0.5,
    }


@dataclass
class ThermalBand:
    """斷面溫度分層：深度範圍（自頂面）、面積、平均溫度。"""
    y_top: float
    y_bot: float
    area: float
    T_mean: float

    @property
    def y_c(self) -> float:
        return (self.y_top + self.y_bot) / 2


@dataclass
class ThermalResult:
    Tu: float                 # 等效均勻溫度 °C
    TL: float                 # 等效線性溫差 °C
    sigma_pos: dict           # 各驗核點 σ_SE（正梯度）MPa，+ 拉
    sigma_neg: dict           # 各驗核點 σ_SE（負梯度 = neg_factor × 正）MPa


def self_equilibrating_stress(bands, I: float, yt: float, h: float, fibers,
                              Ec: float = 30590.0, alpha: float = 1.08e-5,
                              neg_factor: float = -0.30) -> ThermalResult:
    """自平衡溫度應力（斷面積分）。

    bands：ThermalBand list（求 Tu/TL）；I/yt/h：斷面慣性矩/形心距頂/全高；
    fibers：[(名稱, y_自頂, T實際溫度), ...] 驗核點；neg_factor：負梯度係數。
    """
    A_total = sum(b.area for b in bands)
    Tu = sum(b.T_mean * b.area for b in bands) / A_total
    TL = (h / I) * sum(b.T_mean * (yt - b.y_c) * b.area for b in bands)
    Ea = Ec * alpha
    pos, neg = {}, {}
    for name, y, T in fibers:
        T_SE = T - (Tu + TL * (yt - y) / h)
        s = -Ea * T_SE
        pos[name] = s
        neg[name] = neg_factor * s
    return ThermalResult(Tu, TL, pos, neg)


def thermal_service_check(sigma_thermal: float, sigma_base: float,
                          gamma_TG: float = 0.5):
    """Service 疊加：σ_total = σ_base + γ_TG·σ_thermal。

    台灣完全預壓要求 σ ≤ 0（壓）。回傳 (σ_total, ok)。
    """
    total = sigma_base + gamma_TG * sigma_thermal
    return total, total <= 0
