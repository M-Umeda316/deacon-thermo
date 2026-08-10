"""融液の活量モデル。

塩化物融液中の CuCl / CuCl2 の活量を返す。ここが Cu 揮発性の計算の要で、
かつ最も近似の粗い部分。

実装:
  IdealTemkin  : カチオン副格子上の理想混合。錯体形成による安定化を含まないため
                 Cu 蒸気圧を過大評価する。下限見積り用。
  RegularSolution : Temkin + 対相互作用パラメータ。文献の平衡塩素圧データに
                 fit_interactions() でフィットして使う。W は温度依存 (a + bT) も取れる。
                 較正済みのものは calibrated_model()（定数は data.CALIBRATED_INTERACTIONS）。

較正データの向き:
  文献の測定は「仕込み組成（Cu(II)/Cu(I) 比が既知）+ 温度 + p(Cl2)」の形が普通で、
  Cu(II) 比と塩素圧が独立に測られていることは稀。したがって観測は
  ClObservation で保持し、既定は「組成 -> p(Cl2) を予測する」向き（p 型）にしてある。

将来 MQMQA（pycalphad に実装済み）に載せ替える場合は、
ActivityModel のサブクラスを追加するだけでよい。
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq, least_squares

from .data import CALIBRATED_INTERACTIONS, DB
from .species import R

#: 融液成分 -> カチオン
CATIONS = {"CuCl": "Cu+", "CuCl2": "Cu2+", "KCl": "K+"}


def _cation_of(salt: str) -> str:
    return CATIONS.get(salt, salt)


class ActivityModel(ABC):
    @abstractmethod
    def activity(self, salt: str, fractions: dict[str, float], T) -> float:
        ...


class IdealTemkin(ActivityModel):
    """カチオン等価分率をそのまま活量とする。"""

    def activity(self, salt, fractions, T):
        return fractions.get(salt, 0.0)


@dataclass
class RegularSolution(ActivityModel):
    """RT ln(gamma_i) = sum_j W_ij(T) * y_j^2

    W_ij < 0 が安定化（錯体形成）に対応する。

    interactions の値は次のどちらでもよい:

      float          : W = 一定 [J/mol]
      (a, b) タプル  : W(T) = a + b*T [J/mol]（a [J/mol], b [J/mol/K]）

    後者は CALPHAD の L = a + bT と同形で、過剰エントロピー (-b) を分離できる。
    文献の最適化パラメータ（Niazi らの L0 など）をそのまま持ち込むために要る。
    """

    interactions: dict[tuple[str, str], float | tuple[float, float]] = field(
        default_factory=dict
    )

    def _W(self, a: str, b: str, T) -> float:
        w = self.interactions.get((a, b), self.interactions.get((b, a), 0.0))
        if isinstance(w, tuple | list):
            w0, w1 = w
            return w0 + w1 * T
        return w

    def activity(self, salt, fractions, T):
        y_i = fractions.get(salt, 0.0)
        if y_i == 0.0:
            return 0.0
        ln_gamma = sum(
            self._W(salt, other, T) * y_j**2
            for other, y_j in fractions.items()
            if other != salt
        ) / (R * T)
        return y_i * np.exp(ln_gamma)


def calibrated_model() -> RegularSolution:
    """data.CALIBRATED_INTERACTIONS から較正済みモデルを組む。

    既定を差し替えないのは、較正が Cu-K 系（KCl 30 mol%）でしか効いておらず、
    W(Cu 塩化物, LnCl3) が未較正のままだから。Ln を含む系に使うときは
    その未較正分を感度帯として別に振ること（examples/03 参照）。
    """
    return RegularSolution(dict(CALIBRATED_INTERACTIONS))


@dataclass
class Melt:
    """塩化物融液。

    composition : {'CuCl': mol, 'CuCl2': mol, 'KCl': mol, 'SmCl3': mol, ...}
    """

    composition: dict[str, float]
    model: ActivityModel = field(default_factory=IdealTemkin)

    def __post_init__(self):
        self.composition = {k: float(v) for k, v in self.composition.items() if v > 0}

    def cation_fractions(self) -> dict[str, float]:
        total = sum(self.composition.values())
        return {k: v / total for k, v in self.composition.items()}

    def activity(self, salt: str, T) -> float:
        return self.model.activity(salt, self.cation_fractions(), T)


# ---------------------------------------------------------------------------
# Cu(II)/Cu(I) 比を気相の p(Cl2) から決める
# ---------------------------------------------------------------------------

def redox_K(T):
    """2 CuCl(l) + Cl2 = 2 CuCl2(l) の平衡定数。

    融液成分の標準状態は過冷却液体に揃える。redox_split と equilibrium_p_cl2 が
    同じ基準を使うのは、両者が互いの逆関数でなければならないため。
    """
    dG = (
        2 * DB["CuCl2(s)"].G_supercooled_liquid(T)
        - 2 * DB["CuCl(s)"].G_supercooled_liquid(T)
        - DB.G("Cl2(g)", T)
    )
    return float(np.exp(-dG / (R * T)))


def _cu_melt(cu_total, f, diluents=None, model=None) -> Melt:
    """Cu(II) 割合 f の融液を組む。"""
    return Melt(
        {"CuCl": cu_total * (1 - f), "CuCl2": cu_total * f, **dict(diluents or {})},
        model or IdealTemkin(),
    )


def redox_split(cu_total, p_Cl2, T, diluents=None, model=None):
    """2 CuCl(l) + Cl2 = 2 CuCl2(l) を解いて融液を返す。

    Parameters
    ----------
    diluents : {'KCl': mol, 'SmCl3': mol, ...}
    Returns
    -------
    (Cu(II) の割合, Melt)
    """
    diluents = dict(diluents or {})
    model = model or IdealTemkin()
    K = redox_K(T)

    def residual(f):
        m = _cu_melt(cu_total, f, diluents, model)
        a1, a2 = m.activity("CuCl", T), m.activity("CuCl2", T)
        return np.log(a2**2 / (a1**2 * p_Cl2)) - np.log(K)

    f = brentq(residual, 1e-10, 1 - 1e-10)
    return f, _cu_melt(cu_total, f, diluents, model)


def equilibrium_p_cl2(melt: Melt, T) -> float:
    """融液と平衡する Cl2 分圧 [atm]。redox_split の逆向き。

    2 CuCl(l) + Cl2 = 2 CuCl2(l) より p = a(CuCl2)^2 / (a(CuCl)^2 K)。
    文献の平衡塩素圧測定は組成を仕込みで決めて p を測る形なので、
    較正はこちら向きに残差を取るほうが素直になる。
    """
    a1 = melt.activity("CuCl", T)
    a2 = melt.activity("CuCl2", T)
    if a1 <= 0.0:
        raise ValueError("CuCl を含まない融液には Cl2 平衡圧が定義できない")
    return float(a2**2 / (a1**2 * redox_K(T)))


# ---------------------------------------------------------------------------
# 平衡塩素圧データへの較正
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ClObservation:
    """平衡塩素圧の実測点 1 つ。

    Attributes
    ----------
    T        : 温度 [K]。観測ごとに持つので複数温度のデータをそのまま混ぜられる
    cu_total : 融液中の Cu 総量 [mol]（スケールは任意、比だけが効く）
    p_Cl2    : 平衡塩素圧 [atm]
    f_CuII   : Cu(II)/(Cu(I)+Cu(II))。未測定の文献値は None のまま保持できるが、
               どちらの向きにも残差が作れないので fit_interactions が弾く
    diluents : {'KCl': mol, 'LaCl3': mol, ...}
    weight   : 二乗和での重み。p 型（ln 単位）と f 型（0-1 の量）は残差の
               スケールが違うので、混ぜるときはここで揃える
    kind     : "p" -> 組成を条件に p(Cl2) を予測し、残差は ln p_calc - ln p_obs。
                      p は桁で動くので対数で取る。文献データの主用途。
               "f" -> p(Cl2) を条件に Cu(II) 比を解き、残差は f_calc - f_obs。
                      f と p が独立に測られている稀な場合用（旧形式の互換）。
    """

    T: float
    cu_total: float
    p_Cl2: float
    f_CuII: float | None = None
    diluents: dict[str, float] = field(default_factory=dict)
    weight: float = 1.0
    kind: str = "p"

    def __post_init__(self):
        if self.kind not in ("p", "f"):
            raise ValueError(f"kind は 'p' か 'f': {self.kind!r}")
        if self.weight < 0.0:
            raise ValueError("weight は非負")
        object.__setattr__(self, "diluents", dict(self.diluents or {}))

    def melt(self, model: ActivityModel | None = None) -> Melt:
        """組成が既知（f_CuII が与えられている）ときの融液。"""
        if self.f_CuII is None:
            raise ValueError("f_CuII が未測定の観測からは融液を組めない")
        return _cu_melt(self.cu_total, self.f_CuII, self.diluents, model)

    def residual(self, model: ActivityModel) -> float:
        """重み付き残差。least_squares が二乗和を取るので sqrt(weight) を掛ける。"""
        if self.kind == "p":
            r = np.log(equilibrium_p_cl2(self.melt(model), self.T)) - np.log(self.p_Cl2)
        else:
            if self.f_CuII is None:
                raise ValueError("f 型の観測には観測値 f_CuII が必要")
            f_calc, _ = redox_split(
                self.cu_total, self.p_Cl2, self.T, self.diluents, model
            )
            r = f_calc - self.f_CuII
        return float(np.sqrt(self.weight) * r)


def _as_observation(obs, T=None) -> ClObservation:
    """旧形式の 4 つ組 (cu_total, diluents, p_Cl2, f_CuII) も受ける。"""
    if isinstance(obs, ClObservation):
        return obs
    cu_total, diluents, p_obs, f_obs = obs
    if T is None:
        raise ValueError("旧形式の観測には共通温度 T が必要")
    return ClObservation(
        T=T, cu_total=cu_total, p_Cl2=p_obs, f_CuII=f_obs,
        diluents=diluents, kind="f",
    )


def fit_interactions(
    observations, T=None, pairs=None, x0=-5000.0, fixed=None, **least_squares_kw
):
    """実測の平衡塩素圧データに正則溶液パラメータをフィットする。

    新形式::

        fit_interactions([ClObservation(...), ...], pairs)

    旧形式（全点が同一温度で、観測が 4 つ組）::

        fit_interactions([(cu_total, diluents, p_Cl2, f_CuII), ...], T, pairs)

    Parameters
    ----------
    observations : ClObservation のリスト（または旧形式の 4 つ組のリスト）
        CuCl2-CuCl-KCl / CuCl2-CuCl-KCl-LaCl3 融液上の平衡塩素圧測定値を入れる。
    pairs : [('CuCl2', 'KCl'), ('CuCl', 'CuCl2'), ...]
    x0 : スカラー（全パラメータ共通）または pairs と同じ長さの初期値 [J/mol]
    fixed : {('CuCl', 'KCl'): -8000.0 または (-27373.0, 3.94), ...}
        フィットせず固定する W。RegularSolution と同じく (a, b) で
        W = a + bT を渡せる（文献の最適化パラメータの持ち込み用）。

    Notes
    -----
    **フィットする** W は温度非依存（float）として扱う。過剰エントロピー項を
    分離できるほどデータが揃っていないため、温度依存を入れても縮退するだけになる。
    温度依存が独立に判っている組は fixed に (a, b) で渡すこと。

    **塩素圧データは a(CuCl2)^2/a(CuCl)^2 にしか感度がない。** 同じ希釈剤 X に
    対する W(CuCl2, X) と W(CuCl, X) は差だけが決まり、個別には決まらない
    （残差が W(CuCl2,X) - W(CuCl,X) にしか依存しないため、pairs に両方入れると
    正規方程式がランク落ちして初期値依存の答えが返る）。片方は混合熱などの
    独立データから fixed で与えるか、pairs には差を担う組だけを入れること。
    """
    # 新形式では第 2 引数が pairs。数値なら旧形式の T とみなす。
    if pairs is None and T is not None and not np.isscalar(T):
        pairs, T = T, None
    if pairs is None:
        raise ValueError("pairs（相互作用させる成分の組）を指定すること")

    pairs = [tuple(p) for p in pairs]
    obs = [_as_observation(o, T) for o in observations]
    for i, o in enumerate(obs):
        if o.f_CuII is None:
            raise ValueError(f"observations[{i}]: f_CuII 未測定の点はフィットに使えない")

    x0 = np.full(len(pairs), float(x0)) if np.isscalar(x0) else np.asarray(x0, float)
    fixed = dict(fixed or {})

    def build(x) -> RegularSolution:
        return RegularSolution({**fixed, **dict(zip(pairs, x, strict=True))})

    def residual(x):
        model = build(x)
        return [o.residual(model) for o in obs]

    sol = least_squares(residual, x0=x0, **least_squares_kw)
    return build(sol.x), sol
