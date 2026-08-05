"""DH_OXYCHLORIDE 感度解析のテスト。

CLAUDE.md の「系列内の差は絶対値より信頼できる」を、
DH_OXYCHLORIDE を振っても順位が動かないこととして実際に検証する。
"""

import numpy as np
import pytest

from deacon_thermo import (
    chloride_margin_at,
    dHf_oxychloride_threshold,
    flip_threshold,
    flip_thresholds,
    gas_state,
    hydrolysis_margin,
    sweep_margins,
)
from deacon_thermo.data import DB, DH_OXYCHLORIDE, LANTHANIDES

T_OP = 653.15  # 380 C


@pytest.fixture(scope="module")
def gas():
    return gas_state(T_OP, hcl_o2_ratio=2.0)


# --- chloride_margin_at ----------------------------------------------------

def test_margin_at_default_dh_matches_hydrolysis_margin(gas):
    """dh = DH_OXYCHLORIDE では既存の hydrolysis_margin と一致するはず。"""
    for ln in LANTHANIDES:
        assert np.isclose(
            chloride_margin_at(ln, gas, DH_OXYCHLORIDE),
            float(hydrolysis_margin(ln, gas)),
            rtol=1e-10,
        )


def test_margin_moves_one_to_one_with_dh(gas):
    """推定式では dh が dHf(LnOCl) に 1:1 で入るので、余裕も 1:1 で動く。"""
    for ln in ["La", "Sm", "Er"]:
        d = chloride_margin_at(ln, gas, -20.0) - chloride_margin_at(ln, gas, -50.0)
        assert np.isclose(d, 30.0, rtol=1e-9)


def test_deeper_dh_favours_oxychloride(gas):
    """LnOCl を安定化する（dh をより負にする）ほど塩化物側の余裕は減る。"""
    assert chloride_margin_at("Sm", gas, -60.0) < chloride_margin_at("Sm", gas, -10.0)


# --- sweep_margins ---------------------------------------------------------

def test_sweep_covers_all_elements_with_finite_margins(gas):
    sweep = sweep_margins(gas)
    assert set(sweep.margins) == set(LANTHANIDES)
    assert sweep.dh_values.min() == -60.0 and sweep.dh_values.max() == -10.0
    for m in sweep.margins.values():
        assert m.shape == sweep.dh_values.shape
        assert np.all(np.isfinite(m))


def test_ranking_is_independent_of_dh(gas):
    """DH_OXYCHLORIDE は系列共通の平行移動なので、系列内の順位は不変。

    ここが破れたら sensitivity.py か data.py の推定式の実装ミスを疑う。
    """
    sweep = sweep_margins(gas)
    reference = sweep.ranking_at(0)
    for i in range(len(sweep.dh_values)):
        assert sweep.ranking_at(i) == reference


def test_flipped_reports_sign_change_within_sweep(gas):
    sweep = sweep_margins(gas)
    flips = sweep.flipped()
    for ln, flipped in flips.items():
        m = sweep.margins[ln]
        assert flipped == (m.min() < 0.0 < m.max())


# --- flip_threshold --------------------------------------------------------

def test_flip_threshold_zeroes_the_margin(gas):
    for ln in LANTHANIDES:
        dh_star = flip_threshold(ln, gas)
        assert abs(chloride_margin_at(ln, gas, dh_star)) < 1e-6


def test_flip_threshold_consistent_with_dHf_threshold(gas):
    """stability.dHf_oxychloride_threshold（dHf 版の閾値）と整合すること。

    dHf(LnOCl) = 1/3 dHf(Ln2O3) + 1/3 dHf(LnCl3) + dh なので、
    dh の閾値 = dHf の閾値 - 1/3 dHf(Ln2O3) - 1/3 dHf(LnCl3)。
    """
    for ln in ["La", "Sm", "Er"]:
        expected = (
            float(dHf_oxychloride_threshold(ln, gas))
            - DB[f"{ln}2O3(s)"].dHf298 / 3
            - DB[f"{ln}Cl3(s)"].dHf298 / 3
        )
        assert np.isclose(flip_threshold(ln, gas), expected, atol=1e-6)


def test_sm_verdict_is_fragile_but_la_is_not(gas):
    """現行データのピン留め: Sm の判定は ±30 kJ/mol の不確かさ内で反転しうるが、
    La は掃引範囲 (-10, -60) の外でしか反転しない。
    """
    thresholds = flip_thresholds(gas)
    assert -60.0 < thresholds["Sm"] < -10.0
    assert thresholds["La"] < -60.0
