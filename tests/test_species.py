"""既知の熱力学量でピン留めするテスト。

熱力学コードは静かに間違うので、データを差し替えたときに
他が壊れていないかをここで検出する。
"""

import numpy as np
import pytest

from deacon_thermo import DB
from deacon_thermo.species import T_REF


def test_G_at_298_reduces_to_dHf_minus_TS():
    """298.15 K では G = dHf298 - T*S298 に一致するはず。"""
    for name in ["HCl(g)", "H2O(g)", "CuCl(s)", "SmCl3(s)"]:
        sp = DB[name]
        expected = sp.dHf298 * 1000 - T_REF * sp.S298
        assert np.isclose(sp.G(T_REF), expected, rtol=1e-10)


def test_elements_have_zero_formation_enthalpy():
    for name in ["O2(g)", "Cl2(g)"]:
        assert DB[name].dHf298 == 0.0


def test_hcl_oxidation_is_thermodynamically_favourable():
    """4HCl + O2 = 2Cl2 + 2H2O は 380 C で自発的（K > 1）であること。"""
    from deacon_thermo import equilibrium_constant

    assert equilibrium_constant(653.15) > 1.0
    assert np.isclose(DB["H2O(g)"].dHf298, -241.83, atol=0.1)
    assert np.isclose(DB["HCl(g)"].dHf298, -92.31, atol=0.1)


def test_cp_integration_is_consistent():
    """H と S の Cp 積分が整合していること（dH/dT = Cp, T dS/dT = Cp）。"""
    sp = DB["HCl(g)"]
    T, dT = 600.0, 1e-3
    dHdT = (sp.H(T + dT) - sp.H(T - dT)) / (2 * dT)
    TdSdT = T * (sp.S(T + dT) - sp.S(T - dT)) / (2 * dT)
    assert np.isclose(dHdT, TdSdT, rtol=1e-6)


def test_transition_adds_entropy_consistently():
    """相転移で H は dH、S は dH/T_tr だけ跳ぶこと。"""
    sp = DB["CuCl(s)"]
    T_tr, dH_tr = sp.transitions[0]
    eps = 1e-6
    assert np.isclose(sp.H(T_tr + eps) - sp.H(T_tr - eps), dH_tr * 1000, rtol=1e-4)
    assert np.isclose(
        sp.S(T_tr + eps) - sp.S(T_tr - eps), dH_tr * 1000 / T_tr, rtol=1e-4
    )


def test_supercooled_liquid_continuous_at_melting_point():
    """融点で G_liquid と G_solid が一致すること。"""
    sp = DB["CuCl(s)"]
    T_fus = sp.transitions[0][0]
    G_solid = sp.dHf298 * 1000 + sp._cp_int(T_fus) - T_fus * (
        sp.S298 + sp._cp_over_T_int(T_fus)
    )
    assert np.isclose(float(sp.G_supercooled_liquid(T_fus)), G_solid, rtol=1e-9)


def test_supercooled_liquid_above_solid_below_melting_point():
    """融点以下では液体のほうが G が高いこと。"""
    sp = DB["CuCl(s)"]
    T = 653.15
    assert float(sp.G_supercooled_liquid(T)) > float(sp.G(T))


@pytest.mark.parametrize("ln", ["La", "Ce", "Nd", "Sm", "Gd", "Y"])
def test_lanthanide_species_registered(ln):
    for suffix in ["Cl3(s)", "OCl(s)", "2O3(s)"]:
        assert f"{ln}{suffix}" in DB


def test_k3smcl6_is_a_high_temperature_phase():
    """K3SmCl6 の安定下限が 606 K(Seifert Table 8 から導出)を跨ぐこと。

    573 K では KCl + K2SmCl5 に分解し、653 K(380 C の操業点)では安定。
    参照系の操業温度が K3SmCl6 の安定域内にあることは複塩仮説の要。
    """
    for T, stable in [(573.15, False), (653.15, True)]:
        dG = DB.G("K3SmCl6(s)", T) - DB.G("KCl(s)", T) - DB.G("K2SmCl5(s)", T)
        assert (float(dG) < 0) == stable


def test_double_salt_formation_reproduces_seifert_dG573():
    """登録した複塩の生成 dG(573 K) が Seifert Table 8 の dG573 列を再現すること。"""
    for name, n_k, n_l, ln, dg573 in [
        ("K2SmCl5(s)", 2, 1, "Sm", -47.8),
        ("K2GdCl5(s)", 2, 1, "Gd", -51.8),
        ("KSm2Cl7(s)", 1, 2, "Sm", -24.6),  # 表の K0.5SmCl3.5 行 (-12.3) の 2 倍
    ]:
        T = 573.0
        dG = (
            DB.G(name, T) - n_k * DB.G("KCl(s)", T) - n_l * DB.G(f"{ln}Cl3(s)", T)
        ) / 1000
        assert np.isclose(float(dG), dg573, atol=0.15)


def test_oxychloride_estimator_is_consistent():
    """LnOCl の G が 1/3 Ln2O3 + 1/3 LnCl3 + (DH_i, DS_i) の定義どおりであること。

    元素別値(LNOCL_PARAMS)があればそれを、無ければ系列共通既定値を使う。
    """
    from deacon_thermo.data import DH_OXYCHLORIDE, DS_OXYCHLORIDE, LNOCL_PARAMS

    for ln in ["La", "Sm", "Gd", "Dy"]:
        dh, ds = LNOCL_PARAMS.get(ln, (DH_OXYCHLORIDE, DS_OXYCHLORIDE))[:2]
        ds = DS_OXYCHLORIDE if ds is None else ds
        sp = DB[f"{ln}OCl(s)"]
        assert np.isclose(
            sp.dHf298, DB[f"{ln}2O3(s)"].dHf298 / 3 + DB[f"{ln}Cl3(s)"].dHf298 / 3 + dh
        )
        assert np.isclose(
            sp.S298, DB[f"{ln}2O3(s)"].S298 / 3 + DB[f"{ln}Cl3(s)"].S298 / 3 + ds
        )


def test_estimates_are_flagged():
    """推定値が確実に ESTIMATE として印されていること。

    これが緩むと、検証していない数字が黙って結論に入り込む。
    """
    flagged = {s.name for s in DB.needs_verification()}
    assert "CuCl2(g)" in flagged  # JANAF に無く、現在の支配蒸気種。最重要の要検証
    assert "Cu2OCl2(s)" in flagged
    assert "DyOCl(s)" in flagged  # 元素別実測が無く系列共通推定のまま
    # 実測で確定済みの種。EST に戻っていたらデータ退行を疑う
    assert "Cu3Cl3(g)" not in flagged  # JANAF Cl-132
    assert "SmOCl(s)" not in flagged  # Koch & Cunningham 1953 の平衡実測
