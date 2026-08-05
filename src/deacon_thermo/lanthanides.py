"""ランタノイド系列の比較。

設計軸は二つあり、混ぜると解釈できなくなる:

  軸A: Ln 自身が酸化還元に参加するか
       Ce, Pr, Tb は +4、Eu は +2 に届くため第二の酸化還元中心になりうる。
       La, Nd, Gd, Dy, Ho, Er, Y は不活性で、融液化学とルイス酸性のみで効く。
  軸B: イオン半径（ランタノイド収縮）
       塩化物イオン活量、錯体安定性、液相線、Cu 活量係数を連続的に動かす。

このモジュールは軸Bを定量化し、軸Aはフラグとして持つ。

重要な前提
----------
LnOCl の生成エンタルピーは推定値だが、推定誤差は Ln 間で強く相関するため
**系列内の差**は絶対値よりずっと信頼できる。よって
`hydrolysis_ranking()` の順序は使ってよいが、`stable_chloride()` の
絶対判定は単独では信用しないこと。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .data import (
    ACCESSIBLE_OXIDATION_STATES,
    IONIC_RADIUS,
    LANTHANIDES,
    LNCL3_MELTING_POINT,
)
from .gas import GasState
from .stability import dHf_oxychloride_threshold, hydrolysis_dG, hydrolysis_margin


@dataclass(frozen=True)
class LanthanideDescriptors:
    """小データ回帰用の記述子一式。n≈15 なので GNN ではなくこれで足りる。"""

    element: str
    ionic_radius: float  # Å, Shannon CN=6
    redox_active: bool  # 軸A
    accessible_states: tuple[int, ...]
    hydrolysis_dG: float  # kJ/mol, LnCl3 + H2O = LnOCl + 2HCl
    chloride_margin: float  # kJ/mol, 正なら LnCl3 が安定
    dHf_threshold: float  # kJ/mol, LnOCl がこれより負なら反転
    melting_point: float  # K, LnCl3

    def as_row(self) -> dict:
        return {
            "element": self.element,
            "r_ionic": self.ionic_radius,
            "redox_active": int(self.redox_active),
            "dG_hydrolysis": self.hydrolysis_dG,
            "chloride_margin": self.chloride_margin,
            "dHf_threshold": self.dHf_threshold,
            "T_melt_LnCl3": self.melting_point,
        }


def descriptors(ln: str, gas: GasState) -> LanthanideDescriptors:
    states = ACCESSIBLE_OXIDATION_STATES.get(ln, ())
    # Sm(II) は酸化雰囲気では実質届かないので軸Aでは不活性扱いにする
    effective = tuple(s for s in states if not (ln == "Sm" and s == 2))
    return LanthanideDescriptors(
        element=ln,
        ionic_radius=IONIC_RADIUS[ln],
        redox_active=bool(effective),
        accessible_states=states,
        hydrolysis_dG=float(hydrolysis_dG(ln, gas.T)) / 1000,
        chloride_margin=float(hydrolysis_margin(ln, gas)),
        dHf_threshold=float(dHf_oxychloride_threshold(ln, gas)),
        melting_point=float(LNCL3_MELTING_POINT.get(ln, np.nan)),
    )


def survey(gas: GasState, elements=LANTHANIDES) -> list[LanthanideDescriptors]:
    return [descriptors(ln, gas) for ln in elements]


def hydrolysis_ranking(gas: GasState, elements=LANTHANIDES):
    """塩化物として残りやすい順に並べる。系列内の差なので比較的頑健。"""
    rows = survey(gas, elements)
    return sorted(rows, key=lambda d: -d.chloride_margin)


def radius_controls(elements=LANTHANIDES, tol=0.01) -> list[tuple[str, str]]:
    """イオン半径がほぼ同じ Ln のペアを返す。

    Y(0.900) と Ho(0.901) が代表例。4f を持たない Y と 4f を持つ Ho の
    性能が一致すれば効果は純粋にイオン半径・ルイス酸性、
    違えば 4f が関与している、という切り分けができる。
    """
    pairs = []
    els = list(elements)
    for i, a in enumerate(els):
        for b in els[i + 1:]:
            if abs(IONIC_RADIUS[a] - IONIC_RADIUS[b]) <= tol:
                pairs.append((a, b))
    return pairs


def to_dataframe(gas: GasState, elements=LANTHANIDES):
    """pandas がある環境向け。無ければ dict のリストを使うこと。"""
    import pandas as pd

    return pd.DataFrame([d.as_row() for d in survey(gas, elements)])
