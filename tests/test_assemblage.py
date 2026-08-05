"""多金属の安定相集合（カチオン組成断面）のテスト。"""

import numpy as np
import pytest

from deacon_thermo import copper_phases, gas_state, hydrolysis_margin, stability_map
from deacon_thermo.assemblage import (
    Assemblage,
    assemblage_at,
    cation_content,
    cation_grid,
    default_candidates,
    stable_assemblage,
)

T_OP = 653.15  # 380 C


@pytest.fixture(scope="module")
def gas():
    return gas_state(T_OP, hcl_o2_ratio=2.0)


def log_p(gas):
    return float(np.log10(gas.p_O2)), float(np.log10(gas.p_Cl2))


# --- 単一金属の端点が既存実装と一致 ----------------------------------------

def test_single_copper_matches_stability_map(gas):
    """Cu だけの断面は Kellogg 図の操業点と同じ相を返すはず。

    stability_map は格子上の argmin なので phase_at は最近傍に丸まるが、
    操業点では次点（Cu2OCl2）と 8 kJ/mol 離れており丸め誤差では動かない。
    """
    smap = stability_map(copper_phases(), T_OP)
    expected = smap.phase_at(*log_p(gas))

    asm = stable_assemblage({"Cu": 1.0}, gas)
    assert asm.label == (f"{expected}(s)",)


@pytest.mark.parametrize("ln", ["La", "Er"])
def test_lanthanide_matches_hydrolysis_margin(ln, gas):
    """grand potential 経由の判定と加水分解経由の判定が一致すること。

    hydrolysis_margin は p(HCl)^2/p(H2O) で測り、こちらは mu_Cl2/mu_O2 で
    測る。気相が Deacon 平衡 (4HCl + O2 = 2Cl2 + 2H2O) にあれば
    2*mu_HCl - mu_H2O = 2*mu_Cl - mu_O なので両者は厳密に同じ量になる。
    gas フィクスチャは平衡状態なのでこの前提が成り立つ。
    """
    asm = stable_assemblage({ln: 1.0}, gas)
    margin = hydrolysis_margin(ln, gas)
    expected = f"{ln}Cl3(s)" if margin > 0 else f"{ln}OCl(s)"
    assert asm.label == (expected,)
    # La は塩化物側、Er はオキシ塩化物側、という既知の系列傾向も固定しておく
    assert (ln == "La") == (margin > 0)


def test_grand_potential_difference_equals_hydrolysis_margin(gas):
    """二つの経路が一致することを数値そのもので確認する（符号だけでなく値）。"""
    lo2, lcl2 = log_p(gas)
    for ln in ["La", "Sm", "Er"]:
        chloride = assemblage_at(
            {ln: 1.0}, T_OP, lo2, lcl2, [f"{ln}Cl3(s)"]
        ).omega
        oxychloride = assemblage_at(
            {ln: 1.0}, T_OP, lo2, lcl2, [f"{ln}OCl(s)"]
        ).omega
        assert np.isclose(
            (oxychloride - chloride) / 1000, hydrolysis_margin(ln, gas), rtol=1e-9
        )


# --- LP としての健全性 ------------------------------------------------------

def test_cation_balance_is_closed(gas):
    metals = {"Cu": 0.4, "K": 0.35, "Sm": 0.25}
    asm = stable_assemblage(metals, gas)
    balance = asm.metal_balance()
    for el, b in metals.items():
        assert np.isclose(balance[el], b, rtol=1e-9, atol=1e-12)


def test_amounts_are_non_negative(gas):
    asm = stable_assemblage({"Cu": 0.4, "K": 0.35, "Sm": 0.25}, gas)
    assert asm.phases
    assert all(v > 0 for v in asm.phases.values())


def test_phase_count_at_most_number_of_metals(gas):
    """LP の基底解なので共存相数は金属数（等式制約数）以下になる。"""
    for metals in [
        {"Cu": 1.0, "K": 1.0, "Sm": 1.0},
        {"Cu": 0.9, "K": 0.05, "Sm": 0.05},
        {"Cu": 1.0, "K": 0.0, "Sm": 0.0},
    ]:
        asm = stable_assemblage(metals, gas)
        assert len(asm.phases) <= len(metals)


def test_potassium_always_appears_as_kcl(gas):
    """KCl は K の唯一の候補なので、K がある限り必ず立つ。"""
    asm = stable_assemblage({"Cu": 0.5, "K": 0.3, "Sm": 0.2}, gas)
    assert "KCl(s)" in asm.phases
    assert np.isclose(asm.phases["KCl(s)"], 0.3, rtol=1e-9)


def test_zero_amount_metal_gives_binary_section(gas):
    """b=0 の金属を含んでも解け、その金属の相は出てこないこと。"""
    asm = stable_assemblage({"Cu": 0.6, "K": 0.0, "Sm": 0.4}, gas)
    assert "KCl(s)" not in asm.phases
    assert all("Sm" in name or "Cu" in name for name in asm.phases)


def test_degenerate_potentials_are_reported(gas):
    """縮退（同じ Omega の別解）で落ちず、alternatives に出ること。

    候補を CuCl と CuCl2 に絞り、両者の Omega が一致する塩素圧を解く:
      Omega(CuCl) = Omega(CuCl2)  <=>  mu_Cl2 = 2*(G(CuCl2) - G(CuCl))
    これは CuCl2 = CuCl + 1/2 Cl2 の平衡塩素圧そのもの（触媒サイクルの本体）。
    """
    from deacon_thermo.data import DB
    from deacon_thermo.species import R

    lo2, _ = log_p(gas)
    mu_Cl2 = 2 * (DB.G("CuCl2(s)", T_OP) - DB.G("CuCl(s)", T_OP))
    lcl2_tie = (mu_Cl2 - DB.G("Cl2(g)", T_OP)) / (R * T_OP * np.log(10))

    asm = assemblage_at(
        {"Cu": 1.0}, T_OP, lo2, lcl2_tie, ["CuCl(s)", "CuCl2(s)"], degeneracy_tol=1e-3
    )
    assert len(asm.phases) == 1
    assert asm.alternatives  # 相手側が同じ Omega で控えている
    assert set(asm.phases) | set(asm.alternatives) == {"CuCl(s)", "CuCl2(s)"}


def test_unconstrained_metal_is_excluded(gas):
    """収支を課さない金属を含む相を残すと LP が非有界になるので除外される。"""
    asm = stable_assemblage({"Cu": 1.0}, gas, candidates=default_candidates("Sm"))
    assert all(set(cation_content(name)) == {"Cu"} for name in asm.phases)


def test_missing_candidate_raises(gas):
    with pytest.raises(ValueError):
        stable_assemblage({"Cu": 1.0}, gas, candidates=["KCl(s)"])


# --- カチオン組成断面 ------------------------------------------------------

def test_cation_grid_solves_everywhere(gas):
    grid = cation_grid(gas, "Sm", n=12)
    assert len(grid.assemblages) == 13 * 14 // 2
    assert all(isinstance(a, Assemblage) for a in grid.assemblages)
    assert np.allclose(grid.coords.sum(axis=1), 1.0)
    assert all(a.phases for a in grid.assemblages)


def test_cation_grid_vertices_match_single_metal(gas):
    grid = cation_grid(gas, "Sm", n=10)
    for col, el in enumerate(grid.metals):
        vertex = int(np.argmax(grid.coords[:, col]))
        assert np.isclose(grid.coords[vertex, col], 1.0)
        assert grid.labels[vertex] == stable_assemblage({el: 1.0}, gas).label


def test_cation_grid_is_uniform_without_double_salts(gas):
    """複塩が DB に無いうちは各金属が独立に相を選ぶので、内部は一様になる。

    KCuCl3 や K3LnCl6 を登録したらこのテストは分割を検出して落ちるはずで、
    そのときは「どう分かれたか」を確認して期待値を更新すること。
    """
    grid = cation_grid(gas, "Sm", n=8)
    interior = [
        lab for lab, x in zip(grid.labels, grid.coords, strict=True) if (x > 0).all()
    ]
    assert len(set(interior)) == 1
    assert len(next(iter(set(interior)))) == 3
