"""預鑄節塊工法（H5 SBS/BCM）與節塊接縫設計（H6）：節塊重、拼裝期接縫壓、
剪力鍵承載與驗核、黏結 PT 比例。

對應算例_預鑄節塊工法施工設計（H5）、算例_節塊接縫設計（H6）、對應公式卡。
節段橋特徵：節塊在工廠配對澆置→現場拼裝；接縫面**無握裹鋼筋過縫**，剪力由
剪力鍵（道示：設計承載 V_fuk·ξ₁ξ₂Φ）承擔，並要求接縫全程受壓（拼裝期 ≥0.21 MPa）。
介面工程單位：Ac[m² 或 mm²，標於介面]、力[kN]、應力[MPa]。
"""

# 節段橋規範門檻
JOINT_MIN_COMPRESSION_MPa = 0.21   # 拼裝期接縫最小壓應力（AASHTO/道示 ≈30 psi）
BONDED_PT_MIN_RATIO = 0.30         # 內置黏結 PT 最小比例（外置為主時）


def segment_weight(Ac_m2: float, L_seg_m: float, gamma_kNpm3: float = 25.0) -> float:
    """預鑄節塊自重 W = Ac·L_seg·γ，kN（運輸/起重限制典型 ≤30 t ≈ 294 kN）。"""
    return Ac_m2 * L_seg_m * gamma_kNpm3


def joint_min_prestress(Ac_mm2: float, sigma_min: float = JOINT_MIN_COMPRESSION_MPa) -> float:
    """接縫拼裝期所需最小臨時預力 = σ_min·Ac，kN。"""
    return sigma_min * Ac_mm2 / 1e3


def joint_compression(P_kN: float, Ac_mm2: float) -> float:
    """接縫面均勻壓應力 σ = P/Ac，MPa（置中臨時 PT → 均勻）。"""
    return P_kN * 1e3 / Ac_mm2


def shear_key_design_capacity(V_fuk_kN: float, xi_factor: float) -> float:
    """剪力鍵設計承載力 = V_fuk·ξ₁ξ₂Φ（道示；台形鍵 ξ≈0.439、多段鍵≈0.338），kN/鍵。"""
    return V_fuk_kN * xi_factor


def shear_key_utilization(V_sd_kN: float, n_keys: int, V_key_design_kN: float) -> float:
    """剪力鍵驗核比 = (V_sd / n_keys) / V_key_design（≤1 通過）。"""
    return (V_sd_kN / n_keys) / V_key_design_kN


def bonded_pt_ratio(P_bonded_kN: float, P_total_kN: float) -> float:
    """黏結 PT 比例 = P_bonded / P_total（節段橋要求 ≥ BONDED_PT_MIN_RATIO）。"""
    return P_bonded_kN / P_total_kN
