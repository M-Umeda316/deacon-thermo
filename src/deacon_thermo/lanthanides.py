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

記述子の欠測について
--------------------
`k3lncl6_limit` と `chloride_fraction` は K-Ln 複塩データ（Seifert 2002
Table 8）に依存し、La-Gd しか登録が無い。未登録の Ln では前者が None、
後者が 0.0 になるが、いずれも**「複塩が存在しない」ではなく「データが無い」**
を意味する。回帰に載せるときは欠測として扱うこと。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .assemblage import stable_assemblage
from .data import (
    ACCESSIBLE_OXIDATION_STATES,
    DB,
    DH_OXYCHLORIDE,
    DS_OXYCHLORIDE,
    IONIC_RADIUS,
    K_LN_DOUBLE_SALTS,
    LANTHANIDES,
    LNCL3_MELTING_POINT,
    LNOCL_PARAMS,
)
from .gas import GasState
from .stability import dHf_oxychloride_threshold, hydrolysis_dG, hydrolysis_margin

#: 記述子を評価する既定のカチオン組成（参照系 Cu-K-Sm/γ-Al2O3 に相当）。
#: examples/03_series_prediction.py と同じ値。
REFERENCE_CATIONS: dict[str, float] = {"Cu": 0.4, "K": 0.35, "Ln": 0.25}


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
    dh_oxychloride: float  # kJ/mol, 1/3 Ln2O3 + 1/3 LnCl3 -> LnOCl の dH（元素別）
    ds_oxychloride: float  # J/mol/K, 同反応の dS（元素別）
    lnocl_confidence: str  # 上記 2 つの信頼度。回帰の重み付けに使う
    k3lncl6_limit: float | None  # K, K3LnCl6 の安定下限。複塩未登録なら None
    chloride_fraction: float | None  # 参照組成で Ln が塩化物系にいる分率

    def as_row(self) -> dict:
        return {
            "element": self.element,
            "r_ionic": self.ionic_radius,
            "redox_active": int(self.redox_active),
            "dG_hydrolysis": self.hydrolysis_dG,
            "chloride_margin": self.chloride_margin,
            "dHf_threshold": self.dHf_threshold,
            "T_melt_LnCl3": self.melting_point,
            "dH_LnOCl": self.dh_oxychloride,
            "dS_LnOCl": self.ds_oxychloride,
            "LnOCl_confidence": self.lnocl_confidence,
            "T_K3LnCl6_min": self.k3lncl6_limit,
            "chloride_fraction": self.chloride_fraction,
        }


def oxychloride_params(ln: str) -> tuple[float, float, str]:
    """LnOCl 生成反応の (dH [kJ/mol], dS [J/mol/K], 信頼度)。

    元素別実測が無い Ln は系列共通の既定値に落ち、信頼度は estimate になる。
    dH は元素で -35〜-64 kJ/mol とばらつくので、系列一定を仮定しないこと
    （data.py 冒頭の注記）。
    """
    entry = LNOCL_PARAMS.get(ln)
    if entry is None:
        return DH_OXYCHLORIDE, DS_OXYCHLORIDE, "estimate"
    dh, ds, conf, _ = entry
    return float(dh), float(DS_OXYCHLORIDE if ds is None else ds), str(conf.value)


def k3lncl6_stability_limit(ln: str) -> tuple[float, str] | None:
    """K3LnCl6 の安定下限温度 [K] と、その由来を返す。

    Seifert の表に L/H 転移がある Ln は "L"（低温形が別にあるので「この化合物」
    の下限は転移点）、無ければ合成反応 n KCl + LnCl3 -> K3LnCl6 の dG = 0 を
    解いて "S"（それ以下では二元塩化物に分解）。Cp が Neumann-Kopp なので
    dG は厳密に dH - T dS で、根は解析的に出る。

    複塩が data.py に未登録の Ln は None。**「複塩が存在しない」ではなく
    「データが無い」を意味する**（重希土の K 系複塩は Seifert 2002 Table 8 に
    無く、Y については Seifert & Büchel 1998 が未入手）。

    注意
    ----
    ここで返すのは二元塩化物に対する下限であって、実際に相集合で競合するのは
    K2LnCl5 + KCl のほうである（K3LnCl6 は K を 3 使うので K が足りないと
    そちらに負ける）。その交差温度は Sm 609 K / Nd 714 K で、この関数の値
    （Sm 606 K / Nd 169 K）とは意味も大小関係も違う。相集合を知りたいときは
    `assemblage.stable_assemblage` を直接使うこと。
    """
    entry = K_LN_DOUBLE_SALTS.get(f"K3{ln}Cl6(s)")
    if entry is None:
        return None
    _, _, _, dh, ds, transitions = entry
    if transitions:
        return float(transitions[0][0]), "L"
    if ds == 0.0:
        return None
    return dh * 1000.0 / ds, "S"


def chloride_fraction(
    ln: str, gas: GasState, cations: dict[str, float] = REFERENCE_CATIONS
) -> float | None:
    """与えたカチオン組成で、Ln のうち塩化物系の相にいる分率。

    塩化物系 = LnCl3 と K-Ln 複塩（K2LnCl5 / K3LnCl6 / KLn2Cl7 / K3Ln5Cl18）。
    残りは LnOCl（および Ln2O3）に行く。`assemblage.stable_assemblage` の
    線形計画をそのまま使うので、判定は気相が固定するポテンシャルの下で厳密。

    注意
    ----
    **複塩が未登録の Ln では 0.0 になるが、それは「データが無い」であって
    「複塩が存在しない」ではない。** Dy 以降の重希土は Seifert 2002 Table 8 に
    K 系複塩の記載が無く、Y も未入手文献待ち（CLAUDE.md の作業優先順位 2）。
    したがって 0.0 の値を軸B のトレンドに素直に載せてはいけない。
    """
    if ln not in LANTHANIDES:
        return None
    metals = {"Cu": cations["Cu"], "K": cations["K"], ln: cations["Ln"]}
    asm = stable_assemblage(metals, gas)
    total = chloride = 0.0
    for name, amount in asm.phases.items():
        species = DB[name]
        nu = species.elements.get(ln, 0.0)
        if not nu:
            continue
        total += amount * nu
        if "O" not in species.elements:
            chloride += amount * nu
    return chloride / total if total > 0 else None


def descriptors(
    ln: str, gas: GasState, cations: dict[str, float] = REFERENCE_CATIONS
) -> LanthanideDescriptors:
    states = ACCESSIBLE_OXIDATION_STATES.get(ln, ())
    # Sm(II) は酸化雰囲気では実質届かないので軸Aでは不活性扱いにする
    effective = tuple(s for s in states if not (ln == "Sm" and s == 2))
    dh, ds, conf = oxychloride_params(ln)
    limit = k3lncl6_stability_limit(ln)
    return LanthanideDescriptors(
        element=ln,
        ionic_radius=IONIC_RADIUS[ln],
        redox_active=bool(effective),
        accessible_states=states,
        hydrolysis_dG=float(hydrolysis_dG(ln, gas.T)) / 1000,
        chloride_margin=float(hydrolysis_margin(ln, gas)),
        dHf_threshold=float(dHf_oxychloride_threshold(ln, gas)),
        melting_point=float(LNCL3_MELTING_POINT.get(ln, np.nan)),
        dh_oxychloride=dh,
        ds_oxychloride=ds,
        lnocl_confidence=conf,
        k3lncl6_limit=None if limit is None else limit[0],
        chloride_fraction=chloride_fraction(ln, gas, cations),
    )


def survey(
    gas: GasState, elements=LANTHANIDES, cations: dict[str, float] = REFERENCE_CATIONS
) -> list[LanthanideDescriptors]:
    return [descriptors(ln, gas, cations) for ln in elements]


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


def to_dataframe(
    gas: GasState, elements=LANTHANIDES, cations: dict[str, float] = REFERENCE_CATIONS
):
    """pandas がある環境向け。無ければ dict のリストを使うこと。"""
    import pandas as pd

    return pd.DataFrame([d.as_row() for d in survey(gas, elements, cations)])
