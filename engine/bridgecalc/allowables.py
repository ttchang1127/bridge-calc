"""容許值（值寫死＋規範出處，符合 網頁計算器SOP §八）。應力 MPa，壓為負。

每個值標規範條號，便於送審回溯；改版只改此檔。
"""


def comp_service(fc: float) -> float:
    """使用階段混凝土壓應力上限（AASHTO 5.9.4.2）。"""
    return -0.60 * fc


def comp_service_TW_permanent(fc: float) -> float:
    """台灣 永久載重壓應力上限 0.45 f'c。"""
    return -0.45 * fc


def comp_transfer(fci: float) -> float:
    """台灣 施拉階段壓應力上限 0.55 f'ci。"""
    return -0.55 * fci


def tension_full_prestress() -> float:
    """台灣 完全預壓設計：底緣零拉。"""
    return 0.0


def tension_serviceIII(fc: float) -> float:
    """AASHTO Service III 底緣容許拉應力 ≈ 0.5√f'c（MPa，一般腐蝕）。"""
    return 0.5 * fc ** 0.5
