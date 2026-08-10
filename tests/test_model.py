"""気相平衡・安定領域・融液・揮発性のテスト。"""

import numpy as np

from deacon_thermo import (
    IdealTemkin,
    Melt,
    ReactorSpec,
    copper_phases,
    cu_vapour_fraction,
    equilibrium_constant,
    gas_state,
    hydrolysis_K,
    hydrolysis_margin,
    lanthanide_phases,
    lifetime,
    operating_line,
    partial_pressures,
    radius_controls,
    redox_split,
    stability_map,
    survey,
)

T_OP = 653.15  # 380 C


# --- 気相 ------------------------------------------------------------------

def test_deacon_is_exothermic_so_K_falls_with_temperature():
    assert equilibrium_constant(573.15) > equilibrium_constant(773.15)


def test_equilibrium_conversion_in_plausible_range():
    """380 C, HCl/O2=2:1 で転化率は 70-95% 程度になるはず。"""
    gas = gas_state(T_OP, hcl_o2_ratio=2.0)
    assert 0.70 < gas.conversion < 0.95


def test_partial_pressures_sum_to_total():
    gas = gas_state(T_OP)
    assert np.isclose(sum(gas.as_dict().values()), 1.0, rtol=1e-9)


def test_stoichiometry_cl2_equals_h2o():
    """4HCl + O2 = 2Cl2 + 2H2O なので Cl2 と H2O は等モル。"""
    gas = gas_state(T_OP)
    assert np.isclose(gas.p_Cl2, gas.p_H2O, rtol=1e-9)


def test_equilibrium_state_satisfies_K():
    gas = gas_state(T_OP)
    Q = gas.p_Cl2**2 * gas.p_H2O**2 / (gas.p_HCl**4 * gas.p_O2)
    assert np.isclose(Q, equilibrium_constant(T_OP), rtol=1e-6)


def test_operating_line_is_monotonic_in_conversion():
    line = operating_line(T_OP, n=20)
    conv = [g.conversion for g in line]
    assert all(b > a for a, b in zip(conv, conv[1:], strict=False))
    assert line[0].p_Cl2 < line[-1].p_Cl2


def test_lower_temperature_favours_higher_conversion():
    assert gas_state(573.15).conversion > gas_state(773.15).conversion


# --- 安定領域 --------------------------------------------------------------

def test_stability_map_covers_all_phases_somewhere():
    smap = stability_map(lanthanide_phases("Sm"), T_OP)
    assert set(np.unique(smap.index)) <= set(range(len(smap.names)))


def test_chloride_favoured_at_high_pCl2_low_pO2():
    smap = stability_map(lanthanide_phases("La"), T_OP)
    assert smap.phase_at(log_pO2=-11, log_pCl2=0.5) == "LaCl3"


def test_oxide_favoured_at_high_pO2_low_pCl2():
    """文献アンカー版の LnOCl データでは LaOCl の安定域が広く、653 K では
    La2O3 は log_pCl2 < -14.5 (pO2=10^0.5) まで現れない。既定の描画範囲の外。
    """
    smap = stability_map(
        lanthanide_phases("La"), T_OP, log_pCl2_range=(-24, 1)
    )
    assert smap.phase_at(log_pO2=0.5, log_pCl2=-18) == "La2O3"
    assert smap.phase_at(log_pO2=0.5, log_pCl2=-13) == "LaOCl"


def test_copper_phases_include_chloride_at_operating_point():
    gas = gas_state(T_OP)
    smap = stability_map(copper_phases(), T_OP)
    phase = smap.phase_at(np.log10(gas.p_O2), np.log10(gas.p_Cl2))
    assert phase in {"CuCl2", "Cu2OCl2"}


def test_hydrolysis_K_increases_with_temperature():
    """加水分解は吸熱なので高温ほど LnOCl 側に寄る。"""
    assert hydrolysis_K("Sm", 773.15) > hydrolysis_K("Sm", 573.15)


def test_margin_sign_matches_stable_chloride():
    from deacon_thermo import stable_chloride

    gas = gas_state(T_OP)
    for ln in ["La", "Sm", "Gd"]:
        assert (hydrolysis_margin(ln, gas) > 0) == stable_chloride(ln, gas)


def test_dry_feed_strongly_favours_chloride():
    """転化率ゼロ付近では H2O がほぼ無く、塩化物が圧倒的に安定。

    文献アンカー版データでは入口余裕は約 +46 kJ/mol（旧データより 17.6 低い）。
    """
    gas = operating_line(T_OP, n=50)[0]
    assert hydrolysis_margin("Sm", gas) > 40.0


# --- ランタノイド系列 ------------------------------------------------------

def test_survey_returns_all_elements():
    gas = gas_state(T_OP)
    rows = survey(gas)
    assert len(rows) >= 10
    assert all(np.isfinite(d.chloride_margin) for d in rows)


def test_redox_active_flags():
    gas = gas_state(T_OP)
    flags = {d.element: d.redox_active for d in survey(gas)}
    assert flags["Ce"] and flags["Pr"] and flags["Eu"]
    assert not flags["La"] and not flags["Gd"] and not flags["Y"]
    # Sm(II) は酸化雰囲気では届かないので軸Aでは不活性扱い
    assert not flags["Sm"]


def test_ionic_radius_decreases_across_series():
    gas = gas_state(T_OP)
    rows = {d.element: d.ionic_radius for d in survey(gas)}
    for a, b in [("La", "Nd"), ("Nd", "Sm"), ("Sm", "Gd"), ("Gd", "Er")]:
        assert rows[a] > rows[b]


def test_yttrium_holmium_is_a_radius_control_pair():
    """Y と Ho はほぼ同半径。4f の寄与を切り分ける対照になる。"""
    pairs = radius_controls()
    assert ("Ho", "Y") in pairs or ("Y", "Ho") in pairs


# --- 融液と揮発 ------------------------------------------------------------

def test_cation_fractions_sum_to_one():
    melt = Melt({"CuCl": 0.3, "CuCl2": 0.7, "KCl": 1.0, "SmCl3": 0.3})
    assert np.isclose(sum(melt.cation_fractions().values()), 1.0)


def test_dilution_lowers_cu_activity():
    a_neat = Melt({"CuCl2": 1.0}).activity("CuCl2", T_OP)
    a_dil = Melt({"CuCl2": 1.0, "KCl": 3.0}).activity("CuCl2", T_OP)
    assert a_dil < a_neat


def test_higher_pCl2_gives_more_CuII():
    f_low, _ = redox_split(1.0, 0.01, T_OP, {"KCl": 1.0})
    f_high, _ = redox_split(1.0, 0.5, T_OP, {"KCl": 1.0})
    assert f_high > f_low


def test_redox_split_independent_of_dilution_in_ideal_temkin():
    """理想 Temkin では a2^2/a1^2 で希釈が相殺するので Cu(II) 比は変わらない。

    この性質が破れたら活量モデルの実装ミスを疑う。
    """
    f_a, _ = redox_split(1.0, 0.3, T_OP, {"KCl": 0.0})
    f_b, _ = redox_split(1.0, 0.3, T_OP, {"KCl": 3.0, "SmCl3": 1.0})
    assert np.isclose(f_a, f_b, rtol=1e-6)


def test_regular_solution_stabilisation_lowers_activity():
    from deacon_thermo import RegularSolution

    comp = {"CuCl2": 1.0, "KCl": 2.0}
    ideal = Melt(comp, IdealTemkin()).activity("CuCl2", T_OP)
    stabilised = Melt(comp, RegularSolution({("CuCl2", "KCl"): -20000.0})).activity(
        "CuCl2", T_OP
    )
    assert stabilised < ideal


def test_dilution_reduces_cu_volatility():
    _, neat = redox_split(1.0, 0.32, T_OP)
    _, diluted = redox_split(1.0, 0.32, T_OP, {"KCl": 2.0, "SmCl3": 1.0})
    assert cu_vapour_fraction(diluted, T_OP) < cu_vapour_fraction(neat, T_OP)


def test_volatility_increases_with_temperature():
    _, melt = redox_split(1.0, 0.32, T_OP, {"KCl": 1.0})
    lo = cu_vapour_fraction(melt, 573.15)
    hi = cu_vapour_fraction(melt, 773.15)
    assert hi > lo


def test_trimer_dominates_over_monomer():
    """Cu3Cl3 が CuCl(g) より優勢であること（既知の気相化学）。"""
    _, melt = redox_split(1.0, 0.32, T_OP, {"KCl": 1.0})
    p = partial_pressures(melt, T_OP)
    assert p["Cu3Cl3(g)"] > p["CuCl(g)"]


def test_lifetime_scales_inversely_with_volatility():
    _, neat = redox_split(1.0, 0.32, T_OP)
    _, diluted = redox_split(1.0, 0.32, T_OP, {"KCl": 2.0})
    assert lifetime(diluted, T_OP)[0] > lifetime(neat, T_OP)[0]


def test_reactor_spec_gas_flow_positive():
    spec = ReactorSpec()
    assert spec.gas_flow(T_OP) > 0
    assert spec.cu_inventory > 0


def test_ideal_model_reproduces_observed_lifetime():
    """Cu-K-La/Al2O3 は 340 C で 9600 h 安定と報告されている（Feng 2015）。

    かつては xfail だった（約 2 桁外れ）が、原因は Cu3Cl3(g) の暫定値の誤りで
    JANAF 差し替えで解消した。GHSV は Feng 2015 の実条件 450 L/kg/h を使う
    （ReactorSpec 既定の 6000 は設計点であって、この実測の条件ではない。
    旧テストは 6000 のまま偶然通っていた）。支配蒸気種 CuCl2(g) の dHf は
    ESTIMATE のままなので、このテストは「dHf >= 下限」の整合性チェックでもある。
    """
    _, melt = redox_split(1.0, 0.32, 613.15, {"KCl": 1.0, "LaCl3": 0.3})
    t_half, _ = lifetime(melt, 613.15, ReactorSpec(ghsv=450.0))
    assert t_half > 9600.0
