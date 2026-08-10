"""多金属（Cu-K-Ln）の凝縮相集合を、気相が固定したポテンシャルの下で決める。

stability.py は元素 1 種ずつの安定相しか扱わない。ここではカチオンを複数
同時に与え、「どの相が何 mol 共存するか」を返す。

なぜ三元相図ではないか
----------------------
CuCl2-KCl-LnCl3 の三元相図を描こうとすると ill-posed になる（CLAUDE.md）。
実体は Cu-K-Ln-Cl-O-H 系で、Cl と O は気相と自由に交換するため独立変数が
2 つ増えるからである。well-posed な代替が本モジュールの
**カチオン組成断面**で、Cl/O は組成軸ではなくポテンシャルとして扱う:

  固定した (T, mu_Cl2, mu_O2) のもとで、カチオン量 b_M を与えたときの
  grand potential 最小の相集合。

定式化
------
候補相 k（化学量論固体）の式量あたりの量を n_k >= 0 として

    minimize    sum_k n_k * Omega_k
    subject to  sum_k n_k * nu_Mk = b_M      (各金属 M)
                n_k >= 0

    Omega_k = G_k - n_Cl,k * mu_Cl - n_O,k * mu_O      （stability.grand_potential）
    nu_Mk   = 相 k の式量あたりの金属 M の原子数

Cl と O に収支式を課さないのが要点で、気相が無限のリザーバとして働くことを
そのまま表している。H を含む固体は現状 DB に無い（あれば mu_H を足して同型に
扱える）。目的関数も制約も n について線形なので線形計画で厳密に解ける。
汎用の Gibbs 最小化ソルバは不要（CLAUDE.md の方針）。

双対と相律
----------
等式制約の双対変数 y_M が、その相集合における金属 M の化学ポテンシャルに
なる。相 k の被約費用 Omega_k - sum_M y_M * nu_Mk は「その相を析出させる
駆動力」で、解に入っていないのにゼロなら縮退（同じ Omega の別解が存在）。
`Assemblage.alternatives` がこれを拾う。
LP の基底解の性質から、共存相数は金属数以下になる。

複塩について
------------
現状の DB には KCuCl3 や K3LnCl6 といった複塩が無い。複塩が無いうちは
各金属が独立に最安定相を選ぶだけなので、相集合はカチオン比に依存しない。
複塩を data.py に登録して `DOUBLE_SALT_TEMPLATES` か `default_candidates()`
の戻り値に名前を足せば、そのまま候補に入り組成依存が現れる。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog

from .data import DB, LANTHANIDES
from .gas import GasState
from .stability import grand_potential

#: 気相と交換する元素。カチオン収支の対象から外す
EXCHANGED_ELEMENTS = frozenset({"Cl", "O", "H"})

#: Cu の候補相（式量単位。stability.copper_phases() は Cu 1 mol 規格で別物）
COPPER_CANDIDATES = ("CuCl(s)", "CuCl2(s)", "CuO(s)", "Cu2O(s)", "Cu2OCl2(s)")

#: アルカリの候補相。K の候補は KCl のみなので K は必ず KCl として現れる
ALKALI_CANDIDATES = ("KCl(s)",)

#: Ln の候補相テンプレート
LANTHANIDE_TEMPLATES = ("{ln}Cl3(s)", "{ln}OCl(s)", "{ln}2O3(s)")

#: 複塩のテンプレート。DB に登録済みのものだけが候補に入る。
DOUBLE_SALT_TEMPLATES: tuple[str, ...] = (
    "K2{ln}Cl5(s)",
    "K3{ln}Cl6(s)",
    "K{ln}2Cl7(s)",
    "K3{ln}5Cl18(s)",
    "K2CuCl3(s)",
)


def cation_content(name: str) -> dict[str, float]:
    """相 name の式量あたりの金属組成。Cl/O/H は気相と交換するので除く。"""
    return {
        el: n for el, n in DB[name].elements.items() if el not in EXCHANGED_ELEMENTS
    }


def default_candidates(ln: str | None = None) -> list[str]:
    """Cu-K-Ln 系の既定の候補相。

    ln を省略すると Cu-K 系のみ。複塩は DB に登録されているものだけを拾う
    （テンプレートを一括で書いても、一部の Ln しかデータが無い場合に落ちない）。
    """
    names = [*COPPER_CANDIDATES, *ALKALI_CANDIDATES]
    if ln is not None:
        names += [t.format(ln=ln) for t in LANTHANIDE_TEMPLATES]
    for template in DOUBLE_SALT_TEMPLATES:
        if ln is None and "{ln}" in template:
            continue
        name = template.format(ln=ln)
        if name in DB:
            names.append(name)
    return names


@dataclass(frozen=True)
class Assemblage:
    """固定ポテンシャル下で共存する相の集合。

    Attributes
    ----------
    phases : 相名 -> 式量基準の mol。量がゼロの相は含まない
    omega  : 全体の grand potential [J]
    metals : 入力のカチオン量 b_M [mol]
    potentials : 金属 M の化学ポテンシャル [J/mol]（等式制約の双対変数）
    alternatives : 被約費用がゼロだが解に入っていない相（縮退の証拠）
    status : ソルバの終了メッセージ
    """

    phases: dict[str, float]
    omega: float
    metals: dict[str, float]
    potentials: dict[str, float]
    alternatives: tuple[str, ...]
    status: str

    @property
    def label(self) -> tuple[str, ...]:
        """相集合の識別子。ソート済みなので dict のキーに使える。"""
        return tuple(sorted(self.phases))

    def metal_balance(self) -> dict[str, float]:
        """返った相から金属量を再計算する。入力 b と一致するはず。"""
        totals = dict.fromkeys(self.metals, 0.0)
        for name, amount in self.phases.items():
            for el, nu in cation_content(name).items():
                totals[el] = totals.get(el, 0.0) + amount * nu
        return totals


def assemblage_at(
    metals: dict[str, float],
    T: float,
    log_pO2: float,
    log_pCl2: float,
    candidates: list[str],
    amount_tol: float = 1e-9,
    degeneracy_tol: float = 1.0,
) -> Assemblage:
    """(T, log_pO2, log_pCl2) を直接与えるローレベル版。

    Parameters
    ----------
    amount_tol : これ以下の相量は数値ノイズとして捨てる [mol]
    degeneracy_tol : 被約費用がこの範囲 [J/mol] の相を alternatives に入れる。
        既定は実質的な同値判定。大きくすれば「あと何 J/mol で析出するか」の
        近傍探索にも使える。
    """
    if any(b < 0 for b in metals.values()):
        raise ValueError(f"カチオン量が負: {metals}")

    elements = list(metals)
    # 収支を課さない金属を含む相を入れると、その相の n が上に非有界になる
    usable = [k for k in candidates if set(cation_content(k)) <= set(elements)]
    if not usable:
        raise ValueError(f"{elements} に使える候補相が無い")

    covered = {el for k in usable for el in cation_content(k)}
    missing = [el for el, b in metals.items() if b > 0 and el not in covered]
    if missing:
        raise ValueError(f"候補相に含まれない金属: {missing}")

    A = np.array(
        [[cation_content(k).get(el, 0.0) for k in usable] for el in elements],
        dtype=float,
    )
    b = np.array([metals[el] for el in elements], dtype=float)
    c = np.array(
        [grand_potential([(k, 1.0)], T, log_pO2, log_pCl2) for k in usable],
        dtype=float,
    )

    res = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method="highs")
    if not res.success:
        raise ValueError(f"LP が解けなかった: {res.message}")

    n = np.maximum(res.x, 0.0)  # highs は縮退時に -1e-16 程度を返すことがある
    phases = {k: float(v) for k, v in zip(usable, n, strict=True) if v > amount_tol}

    y = np.asarray(res.eqlin.marginals, dtype=float)
    reduced = c - y @ A
    alternatives = tuple(
        k
        for k, r in zip(usable, reduced, strict=True)
        if k not in phases and abs(r) <= degeneracy_tol
    )

    return Assemblage(
        phases=phases,
        omega=float(res.fun),
        metals=dict(metals),
        potentials={el: float(v) for el, v in zip(elements, y, strict=True)},
        alternatives=alternatives,
        status=str(res.message),
    )


def stable_assemblage(
    metals: dict[str, float],
    gas: GasState,
    candidates: list[str] | None = None,
    **kwargs,
) -> Assemblage:
    """気相条件 gas の下でカチオン量 metals に対する安定相集合。

    candidates を省略すると metals に現れる Ln から既定候補を組み立てる。
    """
    if candidates is None:
        candidates = default_candidates(lanthanide_in(metals))
    return assemblage_at(
        metals,
        gas.T,
        float(np.log10(gas.p_O2)),
        float(np.log10(gas.p_Cl2)),
        candidates,
        **kwargs,
    )


def lanthanide_in(metals: dict[str, float]) -> str | None:
    """metals に含まれる Ln を返す。無ければ None、2 種以上なら例外。"""
    found = [el for el in metals if el in LANTHANIDES]
    if len(found) > 1:
        raise ValueError(f"Ln が 2 種以上ある: {found}。candidates を明示すること")
    return found[0] if found else None


# ---------------------------------------------------------------------------
# カチオン組成断面
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CationGrid:
    """Cu-K-Ln カチオン単体の上での相集合分布。

    coords は重心座標（各行の和が 1）なので、そのまま三角図に渡せる。
    描画は呼び出し側の責任で、このモジュールは matplotlib に依存しない。
    """

    metals: tuple[str, ...]
    n: int
    coords: np.ndarray  # (npts, 3) カチオン分率
    assemblages: tuple[Assemblage, ...]
    regions: tuple[tuple[str, ...], ...]  # 一意な相集合ラベル
    region_index: np.ndarray  # (npts,) regions への索引

    @property
    def labels(self) -> list[tuple[str, ...]]:
        return [self.regions[i] for i in self.region_index]

    def region_counts(self) -> dict[tuple[str, ...], int]:
        """相集合ごとの格子点数。断面がどれだけ分割されたかの目安。"""
        counts = np.bincount(self.region_index, minlength=len(self.regions))
        return dict(zip(self.regions, (int(c) for c in counts), strict=True))


def cation_grid(
    gas: GasState,
    ln: str,
    n: int = 60,
    candidates: list[str] | None = None,
) -> CationGrid:
    """Cu:K:Ln 組成単体を n 分割した各点で安定相集合を解く。

    格子点は i+j+k = n の整数点（頂点・辺を含む）で、頂点は単一金属、
    辺は 2 金属断面に対応する。
    """
    metals = ("Cu", "K", ln)
    if candidates is None:
        candidates = default_candidates(ln)

    coords, assemblages = [], []
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            x = np.array([i, j, k], dtype=float) / n
            coords.append(x)
            assemblages.append(
                stable_assemblage(dict(zip(metals, x, strict=True)), gas, candidates)
            )

    labels = [a.label for a in assemblages]
    regions = tuple(sorted(set(labels)))
    lookup = {lab: idx for idx, lab in enumerate(regions)}

    return CationGrid(
        metals=metals,
        n=n,
        coords=np.array(coords),
        assemblages=tuple(assemblages),
        regions=regions,
        region_index=np.array([lookup[lab] for lab in labels], dtype=int),
    )
