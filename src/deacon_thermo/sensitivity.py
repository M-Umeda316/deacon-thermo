"""DH_OXYCHLORIDE の感度解析。

LnOCl の生成エンタルピーは単一パラメータの推定式

    dHf(LnOCl) = 1/3 dHf(Ln2O3) + 1/3 dHf(LnCl3) + DH_OXYCHLORIDE

にぶら下がっている（data.py 参照）。このモジュールはその不確かさが
「塩化物かオキシ塩化物か」の結論をどれだけ揺らすかを定量化する。

要点
----
掃引パラメータ dh は「系列共通既定値 DH_OXYCHLORIDE を dh に置き換える
**共通モードシフト**」で、元素別オフセット（data.py の LNOCL_PARAMS）は
保ったまま全 Ln の chloride_margin を同じ量だけ平行移動する。
したがって系列内の順位は共通モード誤差に対しては不変。ただし実測
（Koch, Yang）で判明した**元素別のばらつき（±15 kJ/mol 級）はこの掃引では
表現されない**ので、順位の頑健性の主張は共通モード誤差に限ること。

DB の数値は書き換えない。LnOCl の dHf298 だけ差し替えた Species を
その場で作って計算する。
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .data import DB, DH_OXYCHLORIDE, LANTHANIDES
from .gas import GasState
from .species import R


def chloride_margin_at(ln: str, gas: GasState, dh_oxychloride: float) -> float:
    """DH_OXYCHLORIDE を仮に dh_oxychloride としたときの塩化物側余裕 [kJ/mol]。

    正なら LnCl3 が安定。dh_oxychloride = DH_OXYCHLORIDE のとき
    stability.hydrolysis_margin() と一致する。
    """
    base = DB[f"{ln}OCl(s)"]
    shifted = replace(base, dHf298=base.dHf298 + (dh_oxychloride - DH_OXYCHLORIDE))
    T = gas.T
    dG = (
        shifted.G(T) + 2 * DB.G("HCl(g)", T)
        - DB.G(f"{ln}Cl3(s)", T) - DB.G("H2O(g)", T)
    )
    return float((R * T * np.log(gas.hydrolysis_quotient) + dG) / 1000)


@dataclass(frozen=True)
class MarginSweep:
    """DH_OXYCHLORIDE 掃引の結果。margins は元素ごとの余裕 [kJ/mol]。"""

    dh_values: np.ndarray
    margins: dict[str, np.ndarray]

    def flipped(self) -> dict[str, bool]:
        """掃引範囲内で判定（塩化物/オキシ塩化物）が反転する Ln。"""
        return {ln: bool(m.min() < 0.0 < m.max()) for ln, m in self.margins.items()}

    def ranking_at(self, i: int) -> tuple[str, ...]:
        """格子点 i での、塩化物として残りやすい順の元素列。"""
        return tuple(
            sorted(self.margins, key=lambda ln: -self.margins[ln][i])
        )


def sweep_margins(
    gas: GasState,
    dh_range: tuple[float, float] = (-10.0, -60.0),
    n: int = 26,
    elements=LANTHANIDES,
) -> MarginSweep:
    """DH_OXYCHLORIDE を dh_range で振り、各 Ln の chloride_margin を返す。"""
    dh_values = np.linspace(dh_range[0], dh_range[1], n)
    margins = {
        ln: np.array([chloride_margin_at(ln, gas, dh) for dh in dh_values])
        for ln in elements
    }
    return MarginSweep(dh_values=dh_values, margins=margins)


def flip_threshold(ln: str, gas: GasState) -> float:
    """判定が反転する DH_OXYCHLORIDE の閾値 [kJ/mol]。

    推定値がこれより負なら LnOCl 側、これより浅ければ LnCl3 側。
    dHf298 は G に加法的に入るので余裕は dh に対して厳密に線形であり、
    2 点の割線で閾値が正確に求まる。
    """
    m0 = chloride_margin_at(ln, gas, DH_OXYCHLORIDE)
    m1 = chloride_margin_at(ln, gas, DH_OXYCHLORIDE + 10.0)
    slope = (m1 - m0) / 10.0
    return DH_OXYCHLORIDE - m0 / slope


def flip_thresholds(gas: GasState, elements=LANTHANIDES) -> dict[str, float]:
    """全 Ln の反転閾値。stability.dHf_oxychloride_threshold の DH 版。"""
    return {ln: flip_threshold(ln, gas) for ln in elements}
