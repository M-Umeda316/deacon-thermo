"""温度依存 W と、較正済み融液モデルの性質。

較正の入力（文献の測定表）はリポジトリ外なので、ここで固定するのは
「登録された定数がそのまま使われるか」と「その定数が意味する物理の向き」の
2 点だけ。フィット手順そのものは test_melt_fit.py が合成データで見ている。
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from deacon_thermo import CALIBRATED_INTERACTIONS, Melt, RegularSolution, calibrated_model
from deacon_thermo.species import R

T_OP = 653.15  # 380 C
T_REF_VOL = 613.15  # 340 C、参照系の寿命条件

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "03_series_prediction.py"


# --- 温度依存 W ------------------------------------------------------------

def test_temperature_dependent_W_follows_a_plus_bT():
    """(a, b) 指定の活量係数が exp((a + bT) y^2 / RT) になること。"""
    a, b = -27373.0, 3.94
    model = RegularSolution({("CuCl", "KCl"): (a, b)})
    fractions = {"CuCl": 0.4, "KCl": 0.6}

    for T in (603.15, 703.15):
        expected = 0.4 * np.exp((a + b * T) * 0.6**2 / (R * T))
        assert np.isclose(model.activity("CuCl", fractions, T), expected, rtol=1e-12)


def test_constant_and_zero_slope_forms_agree():
    """float 指定と (a, 0.0) 指定が完全に一致すること（後方互換の担保）。"""
    fractions = {"CuCl": 0.3, "CuCl2": 0.2, "KCl": 0.5}
    flat = RegularSolution({("CuCl2", "KCl"): -25000.0, ("CuCl", "CuCl2"): 12000.0})
    tupled = RegularSolution(
        {("CuCl2", "KCl"): (-25000.0, 0.0), ("CuCl", "CuCl2"): (12000.0, 0.0)}
    )
    for T in (603.15, T_OP, 703.15):
        for salt in ("CuCl", "CuCl2"):
            assert np.isclose(
                flat.activity(salt, fractions, T),
                tupled.activity(salt, fractions, T),
                rtol=1e-14,
            )


def test_slope_changes_the_temperature_trend():
    """b が効いていること（(a, 0) との差が温度で開く）。"""
    fractions = {"CuCl": 0.4, "KCl": 0.6}
    flat = RegularSolution({("CuCl", "KCl"): (-27373.0, 0.0)})
    sloped = RegularSolution({("CuCl", "KCl"): (-27373.0, 3.94)})
    gaps = [
        abs(sloped.activity("CuCl", fractions, T) - flat.activity("CuCl", fractions, T))
        for T in (603.15, 703.15)
    ]
    assert gaps[1] > gaps[0] > 0.0


def test_missing_pair_is_ideal():
    model = RegularSolution({("CuCl", "KCl"): (-27373.0, 3.94)})
    assert model.activity("KCl", {"KCl": 1.0}, T_OP) == pytest.approx(1.0)


# --- 較正済みモデル --------------------------------------------------------

def test_calibrated_model_uses_the_registered_constants():
    model = calibrated_model()
    assert model.interactions == dict(CALIBRATED_INTERACTIONS)
    # 呼び出しごとに独立（例示スクリプトが W(Cu-Ln) を足しても定数を汚さない）
    model.interactions[("CuCl2", "LaCl3")] = 1.0
    assert ("CuCl2", "LaCl3") not in calibrated_model().interactions
    assert ("CuCl2", "LaCl3") not in CALIBRATED_INTERACTIONS


def test_calibrated_pairs_are_the_three_expected_ones():
    """較正で決まったのはこの 3 組だけ。Ln 系は未較正であることをピン留めする。"""
    assert set(CALIBRATED_INTERACTIONS) == {
        ("CuCl", "KCl"), ("CuCl2", "KCl"), ("CuCl", "CuCl2")
    }
    assert isinstance(CALIBRATED_INTERACTIONS[("CuCl", "KCl")], tuple)


def test_kcl_stabilises_cu_I_more_than_cu_II():
    """KCl 希釈で gamma(CuCl) < 1 < gamma(CuCl2) になること。

    これが較正の物理的な中身そのもの（K2CuCl3 が安定複塩であることと整合）で、
    符号が反転したら結論（K 添加は揮発を増やす向き）ごと変わる。
    """
    model = calibrated_model()
    melt = Melt({"CuCl": 0.35, "CuCl2": 0.15, "KCl": 0.5}, model)
    y = melt.cation_fractions()
    for T in (T_REF_VOL, T_OP):
        assert melt.activity("CuCl", T) / y["CuCl"] < 1.0
        assert melt.activity("CuCl2", T) / y["CuCl2"] > 1.0


def test_more_kcl_deepens_the_effect():
    """希釈量に対して単調（正則溶液なので y^2 で効く）。"""
    model = calibrated_model()
    gammas = []
    for kcl in (0.5, 1.0, 2.0):
        melt = Melt({"CuCl": 0.5, "CuCl2": 0.5, "KCl": kcl}, model)
        y = melt.cation_fractions()
        gammas.append(melt.activity("CuCl", T_OP) / y["CuCl"])
    assert gammas[0] > gammas[1] > gammas[2]


def test_calibration_shortens_the_predicted_lifetime():
    """較正すると参照系の寿命は理想 Temkin より短くなる（gamma(CuCl2) > 1 のため）。

    実測の「>= 9600 h 安定」を理想 11,900 h はまたぐが較正後は下回る。この不一致は
    CuCl2(g) が ESTIMATE であること + W(Cu-LaCl3) 未較正の範囲内なので、
    どちらが正しいかはまだ決められない。値が動いたら理由を追うこと。
    """
    from deacon_thermo import lifetime, redox_split

    diluents = {"KCl": 1.0, "LaCl3": 0.3}
    _, ideal = redox_split(1.0, 0.32, T_REF_VOL, diluents)
    _, fitted = redox_split(1.0, 0.32, T_REF_VOL, diluents, calibrated_model())
    t_ideal = lifetime(ideal, T_REF_VOL)[0]
    t_fitted = lifetime(fitted, T_REF_VOL)[0]
    assert t_fitted < t_ideal
    assert 0.3 < t_fitted / t_ideal < 0.9


# --- 系列予測スクリプトのスモーク ------------------------------------------

@pytest.fixture(scope="module")
def series():
    """examples/03 をパス指定で読み込む（examples はパッケージではない）。"""
    spec = importlib.util.spec_from_file_location("series_prediction", EXAMPLE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_half_life_band_is_ordered(series):
    """W(Cu-Ln) が負（安定化）ほど寿命が延びること。"""
    lo, mid, hi = series.half_lives("Sm")
    assert lo > mid > hi > 0.0


def test_half_lives_are_currently_ln_independent(series):
    """LnCl3 が理想希釈剤にすぎない（W(Cu-Ln) 未較正）ことの可視化。

    ここが破れたら W(Cu-LnCl3) の較正が入ったということなので、
    examples/03 の注意書きを更新すること。
    """
    assert series.half_lives("La") == pytest.approx(series.half_lives("Er"))


def test_k3lncl6_lower_bound_kinds(series):
    """L/H 転移がある Ln は転移点、無い Ln は合成 dG=0 の根を返す。"""
    assert series.k3lncl6_lower_bound("Sm") == pytest.approx((606.1, "L"), rel=1e-9)
    T_ce, kind = series.k3lncl6_lower_bound("Ce")
    assert kind == "S"
    assert 250.0 < T_ce < 400.0
    assert series.k3lncl6_lower_bound("La") is None  # 表に無い（存在しない）


def test_assemblage_summary_accounts_for_all_of_the_lanthanide(series):
    from deacon_thermo import gas_state

    gas = gas_state(series.T_GAS, hcl_o2_ratio=2.0)
    summary = series.assemblage_summary("Sm", gas)
    assert "Sm" in summary
    percentages = [
        float(token.rstrip("%")) for token in summary.split() if token.endswith("%")
    ]
    assert sum(percentages) == pytest.approx(100.0, abs=1.0)
