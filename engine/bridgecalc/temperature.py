"""溫度梯度載重定義（T1，部分）。

對應公式卡 T1、算例_溫度梯度應力設計。
本模組只含**載重定義**（AASHTO 正梯度 Zone 值、PC 箱梁負梯度 −0.30×、γ_TG）。

⚠️ **未含自平衡應力檢核**：σ_SE 需對非矩形斷面（頂板/腹板/底板不同寬度）逐層積分求
等效均勻溫度 Tu 與曲率，且 T1 算例用自含斷面（h=2,000、f'c=35，非配置 A）。屬較大工程，
見 算例_溫度梯度應力設計（負梯度底板 +2.00 MPa 拉為控制工況、底板需有效預壓 ≥ 2.5 MPa）。
"""


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
