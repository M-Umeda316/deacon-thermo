"""融液活量モデルの較正のテスト。

既知の W から擬似観測を作り、fit_interactions が同じ W を取り戻せるかを見る
（合成データでの往復）。実データを入れる前に、較正の向き・残差の定義・
何が決まって何が決まらないかをここで固定しておく。
"""

import numpy as np
import pytest

from deacon_thermo import IdealTemkin, Melt, RegularSolution, redox_split
from deacon_thermo.melt import (
    ClObservation,
    equilibrium_p_cl2,
    fit_interactions,
    redox_K,
)
from deacon_thermo.species import R

T_OP = 653.15  # 380 C

#: 合成データの真値。CuCl2-KCl はクロロ銅酸錯体による安定化を想定して負、
#: CuCl-CuCl2 は Cu(I)/Cu(II) 間の非理想性として正に取ってある。
TRUE_W = {("CuCl2", "KCl"): -25000.0, ("CuCl", "CuCl2"): 12000.0}
FIT_PAIRS = [("CuCl2", "KCl"), ("CuCl", "CuCl2")]

#: (T [K], Cu(II) 割合, KCl [mol])。複数温度・複数組成を混ぜる。
CONDITIONS = [
    (T, f, kcl)
    for T in (603.15, 653.15, 703.15)
    for f in (0.25, 0.5, 0.75)
    for kcl in (0.5, 1.0, 2.0)
]


def synthetic(true_w=None, kind="p", ln_noise=0.0, seed=0):
    """既知の W から擬似観測を作る。ln_noise > 0 なら p に対数正規ノイズを乗せる。"""
    model = RegularSolution(dict(true_w or TRUE_W))
    rng = np.random.default_rng(seed)
    obs = []
    for T, f, kcl in CONDITIONS:
        melt = Melt({"CuCl": 1.0 - f, "CuCl2": f, "KCl": kcl}, model)
        p = equilibrium_p_cl2(melt, T)
        if ln_noise:
            p *= float(np.exp(rng.normal(0.0, ln_noise)))
        obs.append(
            ClObservation(
                T=T, cu_total=1.0, p_Cl2=p, f_CuII=f, diluents={"KCl": kcl}, kind=kind
            )
        )
    return obs


def relative_error(fitted: RegularSolution, true_w=None):
    true_w = dict(true_w or TRUE_W)
    return {
        pair: abs(fitted.interactions[pair] - w) / abs(w) for pair, w in true_w.items()
    }


# --- equilibrium_p_cl2 と redox_split の整合 --------------------------------

@pytest.mark.parametrize(
    "model",
    [IdealTemkin(), RegularSolution(TRUE_W)],
    ids=["ideal", "regular"],
)
def test_equilibrium_p_cl2_inverts_redox_split(model):
    """redox_split で得た融液を戻せば元の p(Cl2) が返ること。

    較正の残差はこの逆写像の上に載るので、ここがずれると全部ずれる。
    """
    for T in (603.15, T_OP, 703.15):
        for p in (1e-5, 1e-3, 0.1):
            _, melt = redox_split(1.0, p, T, {"KCl": 1.0}, model)
            assert np.isclose(equilibrium_p_cl2(melt, T), p, rtol=1e-8)


def test_redox_split_inverts_equilibrium_p_cl2():
    """逆向き（組成 -> p -> 組成）でも一致すること。"""
    model = RegularSolution(TRUE_W)
    for f in (0.2, 0.5, 0.9):
        melt = Melt({"CuCl": 1.0 - f, "CuCl2": f, "KCl": 1.5}, model)
        p = equilibrium_p_cl2(melt, T_OP)
        f_back, _ = redox_split(1.0, p, T_OP, {"KCl": 1.5}, model)
        assert np.isclose(f_back, f, rtol=1e-8)


def test_equilibrium_p_cl2_matches_ideal_closed_form():
    """理想 Temkin では p = (y2/y1)^2 / K に落ちる。"""
    melt = Melt({"CuCl": 0.6, "CuCl2": 0.4, "KCl": 2.0}, IdealTemkin())
    expected = (0.4 / 0.6) ** 2 / redox_K(T_OP)
    assert np.isclose(equilibrium_p_cl2(melt, T_OP), expected, rtol=1e-12)


# --- 残差 ------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["p", "f"])
def test_residual_vanishes_at_true_parameters(kind):
    model = RegularSolution(TRUE_W)
    for o in synthetic(kind=kind):
        assert abs(o.residual(model)) < 1e-8


def test_weight_enters_as_sqrt_in_the_residual():
    """least_squares は二乗和を取るので、残差には sqrt(weight) が掛かる。"""
    model = IdealTemkin()
    base = synthetic()[0]
    heavy = ClObservation(
        T=base.T, cu_total=base.cu_total, p_Cl2=base.p_Cl2, f_CuII=base.f_CuII,
        diluents=base.diluents, weight=4.0,
    )
    assert np.isclose(heavy.residual(model), 2.0 * base.residual(model), rtol=1e-12)


# --- 往復（ノイズ無し） ----------------------------------------------------

@pytest.mark.parametrize("kind", ["p", "f"])
def test_roundtrip_recovers_parameters(kind):
    fitted, sol = fit_interactions(synthetic(kind=kind), FIT_PAIRS)
    assert sol.success
    assert max(relative_error(fitted).values()) < 1e-4


def test_roundtrip_with_mixed_observation_types():
    """p 型と f 型を混ぜても同じ答えに行くこと。"""
    mixed = synthetic(kind="p")[::2] + synthetic(kind="f")[1::2]
    fitted, _ = fit_interactions(mixed, FIT_PAIRS)
    assert max(relative_error(fitted).values()) < 1e-4


def test_single_parameter_fit_matches_analytic_solution():
    """1 パラメータなら手で解ける。

    ln p = 2 ln(y2/y1) + 2 (W2K - W1K) yK^2 / RT - ln K なので
    W2K = W1K + (ln p + ln K - 2 ln(y2/y1)) RT / (2 yK^2)。
    """
    w1k = -8000.0
    true_w = {("CuCl2", "KCl"): -25000.0, ("CuCl", "KCl"): w1k}
    f, kcl = 0.4, 1.0
    melt = Melt({"CuCl": 1.0 - f, "CuCl2": f, "KCl": kcl}, RegularSolution(true_w))
    p = equilibrium_p_cl2(melt, T_OP)
    obs = [
        ClObservation(T=T_OP, cu_total=1.0, p_Cl2=p, f_CuII=f, diluents={"KCl": kcl})
    ]

    y = melt.cation_fractions()
    analytic = w1k + (
        (np.log(p) + np.log(redox_K(T_OP)) - 2 * np.log(y["CuCl2"] / y["CuCl"]))
        * R * T_OP / (2 * y["KCl"] ** 2)
    )
    assert np.isclose(analytic, -25000.0, rtol=1e-9)

    fitted, _ = fit_interactions(
        obs, [("CuCl2", "KCl")], fixed={("CuCl", "KCl"): w1k}
    )
    assert np.isclose(fitted.interactions[("CuCl2", "KCl")], analytic, rtol=1e-6)


def test_chloride_pressure_fixes_only_the_difference_of_the_two_cu_pairs():
    """W(CuCl2,X) と W(CuCl,X) は差しか決まらない（ランク落ち）。

    残差が a(CuCl2)^2/a(CuCl)^2 にしか依存しないため。ここが破れたら
    fit_interactions の docstring の注意書きごと見直すこと。
    """
    true_w = {("CuCl2", "KCl"): -25000.0, ("CuCl", "KCl"): -8000.0}
    obs = synthetic(true_w)
    degenerate = [("CuCl2", "KCl"), ("CuCl", "KCl")]

    fits = [
        fit_interactions(obs, degenerate, x0=x)[0].interactions for x in (0.0, -20000.0)
    ]
    diffs = [w[("CuCl2", "KCl")] - w[("CuCl", "KCl")] for w in fits]
    for d in diffs:
        assert np.isclose(d, -17000.0, rtol=1e-6)
    # 差は決まるが個別の値は初期値で動く
    assert abs(fits[0][("CuCl2", "KCl")] - fits[1][("CuCl2", "KCl")]) > 1000.0

    # 片方を独立データで固定すればもう片方は一意に決まる
    fitted, _ = fit_interactions(
        obs, [("CuCl2", "KCl")], fixed={("CuCl", "KCl"): -8000.0}
    )
    assert np.isclose(fitted.interactions[("CuCl2", "KCl")], -25000.0, rtol=1e-5)


# --- 往復（ノイズ有り） ----------------------------------------------------

def test_noisy_data_still_lands_near_the_truth():
    """ln p に 5% のばらつきを入れても数 % で戻ること。"""
    fitted, _ = fit_interactions(synthetic(ln_noise=0.05, seed=3), FIT_PAIRS)
    err = relative_error(fitted)
    assert max(err.values()) < 0.10


def test_noise_free_fit_is_better_than_noisy_fit():
    """ノイズが誤差として現れていること（往復テストが自明でないことの確認）。"""
    clean, _ = fit_interactions(synthetic(), FIT_PAIRS)
    noisy, _ = fit_interactions(synthetic(ln_noise=0.05, seed=3), FIT_PAIRS)
    assert max(relative_error(clean).values()) < max(relative_error(noisy).values())


# --- 旧形式の互換 ----------------------------------------------------------

def test_legacy_tuple_observations_are_still_accepted():
    """(cu_total, diluents, p_Cl2, f_CuII) + 共通 T の旧呼び出しが動くこと。"""
    true_w = {("CuCl2", "KCl"): -25000.0, ("CuCl", "CuCl2"): 12000.0}
    model = RegularSolution(true_w)
    legacy = []
    for f in (0.25, 0.4, 0.6, 0.75):
        for kcl in (0.5, 1.0, 2.0):
            melt = Melt({"CuCl": 1.0 - f, "CuCl2": f, "KCl": kcl}, model)
            legacy.append((1.0, {"KCl": kcl}, equilibrium_p_cl2(melt, T_OP), f))

    fitted, _ = fit_interactions(legacy, T_OP, FIT_PAIRS)
    assert max(relative_error(fitted, true_w).values()) < 1e-3


def test_unmeasured_cu_split_is_rejected():
    """f_CuII が未測定の点は残差が作れないので、黙って無視せず弾くこと。"""
    obs = synthetic()
    obs.append(ClObservation(T=T_OP, cu_total=1.0, p_Cl2=1e-3, diluents={"KCl": 1.0}))
    with pytest.raises(ValueError, match="f_CuII"):
        fit_interactions(obs, FIT_PAIRS)


def test_invalid_kind_is_rejected():
    with pytest.raises(ValueError, match="kind"):
        ClObservation(T=T_OP, cu_total=1.0, p_Cl2=1e-3, f_CuII=0.5, kind="x")
