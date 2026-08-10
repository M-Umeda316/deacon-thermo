"""拡張した Ln 記述子（複塩・元素別 LnOCl 量）のテスト。

記述子は実験計画（docs/experiment_plan.md）と回帰の入力なので、
data.py との対応が切れていないことと、欠測が「0」に化けていないことを
固定しておく。
"""

import pytest

from deacon_thermo import (
    LANTHANIDES,
    chloride_fraction,
    descriptors,
    gas_state,
    k3lncl6_stability_limit,
    oxychloride_params,
    survey,
)
from deacon_thermo.data import DS_OXYCHLORIDE, LNOCL_PARAMS

T_OP = 653.15  # 380 C


@pytest.fixture(scope="module")
def gas():
    return gas_state(T_OP, hcl_o2_ratio=2.0)


@pytest.fixture(scope="module")
def rows(gas):
    return {d.element: d for d in survey(gas)}


# --- 元素別 LnOCl 量 --------------------------------------------------------

def test_oxychloride_descriptors_match_data(rows):
    """dh/ds/信頼度が LNOCL_PARAMS の転記そのものであること。

    ここがずれると「系列一定を仮定しない」という設計（data.py 冒頭）が
    黙って壊れる。
    """
    for ln in LANTHANIDES:
        dh, ds, conf, _ = LNOCL_PARAMS[ln]
        d = rows[ln]
        assert d.dh_oxychloride == pytest.approx(dh)
        assert d.ds_oxychloride == pytest.approx(DS_OXYCHLORIDE if ds is None else ds)
        assert d.lnocl_confidence == conf.value
        assert d.as_row()["dH_LnOCl"] == pytest.approx(dh)


def test_unknown_element_falls_back_to_series_default():
    """元素別値の無い Ln は系列共通既定値に落ち、信頼度が estimate になること。"""
    dh, ds, conf = oxychloride_params("Tm")  # LNOCL_PARAMS に無い
    assert (dh, ds) == (pytest.approx(-58.0), pytest.approx(DS_OXYCHLORIDE))
    assert conf == "estimate"


def test_dh_is_not_constant_across_series(rows):
    """実測 dH は元素で 25 kJ/mol 以上ばらつく（旧設計の系列一定は誤り）。"""
    values = [rows[ln].dh_oxychloride for ln in LANTHANIDES]
    assert max(values) - min(values) > 25.0


# --- K3LnCl6 の安定下限 ----------------------------------------------------

def test_k3lncl6_limit_sm_is_low_temperature_transition(rows):
    """Sm は L/H 転移（606 K）が下限。操業温度 653 K は上側にある。"""
    assert rows["Sm"].k3lncl6_limit == pytest.approx(606.1, abs=1.0)
    assert k3lncl6_stability_limit("Sm")[1] == "L"
    assert rows["Sm"].k3lncl6_limit < T_OP


def test_k3lncl6_limit_nd_is_synthesis_root(rows):
    """Nd は転移が無く合成 dG=0 の根（169 K）。Sm と由来が違う。"""
    assert k3lncl6_stability_limit("Nd")[1] == "S"
    assert rows["Nd"].k3lncl6_limit == pytest.approx(169.2, abs=1.0)


@pytest.mark.parametrize("ln", ["La", "Dy", "Ho", "Er", "Y"])
def test_k3lncl6_limit_is_none_when_unregistered(ln, rows):
    """K3LnCl6 が未登録の Ln は None（0 ではない）。

    La は K3La5Cl18 と K2LaCl5 は持つが K3LaCl6 は存在せず（Seifert の表にも
    無い）、Dy 以降は重希土で Seifert 2002 Table 8 に記載自体が無い。
    どちらも「複塩が無い」ことの証拠ではないので、None のまま欠測として
    扱えることをここで保証する。
    """
    assert k3lncl6_stability_limit(ln) is None
    assert rows[ln].k3lncl6_limit is None
    assert rows[ln].as_row()["T_K3LnCl6_min"] is None


# --- 塩化物分率（相分配）---------------------------------------------------

def test_chloride_fraction_reference_composition(rows):
    """参照組成での Ln 分配。

    La は K3La5Cl18（K を 0.6/Ln しか使わない）が LaOCl に勝つので全量が
    塩化物側。Sm/Gd などは K2LnCl5 が K を使い切った時点で頭打ちになり
    70 %（K:Ln = 1.4 の予算で決まる）。Er は複塩が未登録で ErOCl 一相。
    """
    assert rows["La"].chloride_fraction == pytest.approx(1.0)
    assert 0.6 < rows["Sm"].chloride_fraction < 0.8
    assert rows["Gd"].chloride_fraction == pytest.approx(rows["Sm"].chloride_fraction)
    assert rows["Er"].chloride_fraction == pytest.approx(0.0)


def test_chloride_fraction_rises_with_potassium(gas, rows):
    """K を増やせば Sm はより多く塩化物として係留される（KCl 添加の役割）。"""
    poor = chloride_fraction("Sm", gas, {"Cu": 0.4, "K": 0.10, "Ln": 0.25})
    rich = chloride_fraction("Sm", gas, {"Cu": 0.2, "K": 0.60, "Ln": 0.20})
    assert poor < rows["Sm"].chloride_fraction < rich
    assert rich == pytest.approx(1.0)


# --- 後方互換 --------------------------------------------------------------

def test_existing_api_preserved(gas):
    """既存の呼び出し形と as_row キーが残っていること。"""
    row = descriptors("Sm", gas).as_row()  # cations 省略で既定値に落ちる
    for key in [
        "element", "r_ionic", "redox_active", "dG_hydrolysis",
        "chloride_margin", "dHf_threshold", "T_melt_LnCl3",
    ]:
        assert key in row
