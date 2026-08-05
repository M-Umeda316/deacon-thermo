"""凝縮相の安定領域図（Kellogg 図）と加水分解平衡。

固定した化学ポテンシャル (mu_Cl2, mu_O2) のもとで、対象元素 1 mol あたりの
grand potential が最小の相が安定になる。候補相が有限個なので、
汎用の Gibbs 最小化ソルバは不要で、格子上の argmin で済む。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .data import DB
from .gas import GasState
from .species import R

LN10 = np.log(10)

#: 対象元素 1 mol あたりに規格化した候補相の定義
Phases = dict[str, list[tuple[str, float]]]


def lanthanide_phases(ln: str) -> Phases:
    """Ln 系の候補相。Ln 1 mol あたり。"""
    return {
        f"{ln}Cl3": [(f"{ln}Cl3(s)", 1.0)],
        f"{ln}OCl": [(f"{ln}OCl(s)", 1.0)],
        f"{ln}2O3": [(f"{ln}2O3(s)", 0.5)],
    }


def copper_phases() -> Phases:
    """Cu 系の候補相。Cu 1 mol あたり。"""
    return {
        "CuCl": [("CuCl(s)", 1.0)],
        "CuCl2": [("CuCl2(s)", 1.0)],
        "CuO": [("CuO(s)", 1.0)],
        "Cu2O": [("Cu2O(s)", 0.5)],
        "Cu2OCl2": [("Cu2OCl2(s)", 0.5)],
    }


def grand_potential(phase_def, T, log_pO2, log_pCl2):
    """対象元素 1 mol あたりの grand potential [J/mol]。

        Omega = G_phase - n_Cl * mu_Cl - n_O * mu_O
        mu_Cl = mu(Cl2)/2,  mu_O = mu(O2)/2
    """
    mu_Cl2 = DB.G("Cl2(g)", T) + R * T * LN10 * log_pCl2
    mu_O2 = DB.G("O2(g)", T) + R * T * LN10 * log_pO2

    omega = 0.0
    for name, coeff in phase_def:
        sp = DB[name]
        omega = omega + coeff * sp.G(T)
        omega = omega - coeff * sp.elements.get("Cl", 0) * mu_Cl2 / 2
        omega = omega - coeff * sp.elements.get("O", 0) * mu_O2 / 2
    return omega


@dataclass
class StabilityMap:
    log_pO2: np.ndarray
    log_pCl2: np.ndarray
    index: np.ndarray  # 各格子点での最安定相のインデックス
    names: list[str]

    def phase_at(self, log_pO2: float, log_pCl2: float) -> str:
        i = int(np.abs(self.log_pO2 - log_pO2).argmin())
        j = int(np.abs(self.log_pCl2 - log_pCl2).argmin())
        return self.names[self.index[j, i]]


def stability_map(phases: Phases, T, log_pO2_range=(-12, 1),
                  log_pCl2_range=(-14, 1), n=400) -> StabilityMap:
    x = np.linspace(*log_pO2_range, n)
    y = np.linspace(*log_pCl2_range, n)
    X, Y = np.meshgrid(x, y)
    names = list(phases)
    stack = np.array([grand_potential(phases[k], T, X, Y) for k in names])
    return StabilityMap(x, y, np.argmin(stack, axis=0), names)


# ---------------------------------------------------------------------------
# 加水分解（O2 に依存しないため、より頑健な指標）
# ---------------------------------------------------------------------------

def hydrolysis_dG(ln: str, T):
    """LnCl3 + H2O = LnOCl + 2 HCl の dG [J/mol]。"""
    return (
        DB.G(f"{ln}OCl(s)", T) + 2 * DB.G("HCl(g)", T)
        - DB.G(f"{ln}Cl3(s)", T) - DB.G("H2O(g)", T)
    )


def hydrolysis_K(ln: str, T):
    """K = p(HCl)^2 / p(H2O)。実測の Q がこれより大きければ塩化物側が安定。"""
    return np.exp(-hydrolysis_dG(ln, T) / (R * T))


def stable_chloride(ln: str, gas: GasState) -> bool:
    """与えられた気相条件で LnCl3 が LnOCl より安定か。"""
    return gas.hydrolysis_quotient > hydrolysis_K(ln, gas.T)


def hydrolysis_margin(ln: str, gas: GasState):
    """塩化物側の安定余裕 [kJ/mol]。正なら LnCl3 が安定。

    これが数十 kJ/mol 未満なら、LnOCl のデータ誤差で結論が反転しうる。
    """
    Q = gas.hydrolysis_quotient
    return (R * gas.T * np.log(Q) + hydrolysis_dG(ln, gas.T)) / 1000


def dHf_oxychloride_threshold(ln: str, gas: GasState):
    """LnOCl が安定側に転ぶ dHf298(LnOCl) の閾値 [kJ/mol]。

    推定値が信用できなくても「閾値のどちら側か」で議論できるようにするための逆算。
    """
    sp = DB[f"{ln}OCl(s)"]
    T = gas.T
    G_crit = (
        DB.G(f"{ln}Cl3(s)", T) + DB.G("H2O(g)", T) - 2 * DB.G("HCl(g)", T)
        - R * T * np.log(gas.hydrolysis_quotient)
    )
    return sp.dHf298 + (G_crit - sp.G(T)) / 1000
