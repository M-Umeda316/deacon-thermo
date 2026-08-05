"""気相の Deacon 平衡。

    4 HCl + O2 = 2 Cl2 + 2 H2O

反応条件（温度・供給比・転化率）から p(HCl), p(O2), p(Cl2), p(H2O) を返す。
凝縮相の安定性はこの4つの分圧で決まるので、ここが全ての入口になる。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .data import DB
from .species import R

REACTION = {"HCl(g)": -4, "O2(g)": -1, "Cl2(g)": 2, "H2O(g)": 2}


@dataclass(frozen=True)
class GasState:
    """気相の状態。分圧は atm。"""

    T: float
    p_HCl: float
    p_O2: float
    p_Cl2: float
    p_H2O: float
    conversion: float
    extent: float

    def as_dict(self) -> dict[str, float]:
        return {
            "HCl(g)": self.p_HCl,
            "O2(g)": self.p_O2,
            "Cl2(g)": self.p_Cl2,
            "H2O(g)": self.p_H2O,
        }

    @property
    def hydrolysis_quotient(self) -> float:
        """Q = p(HCl)^2 / p(H2O)。LnCl3 -> LnOCl の加水分解の反応比。"""
        return self.p_HCl**2 / self.p_H2O


def equilibrium_constant(T):
    """Deacon 反応の平衡定数。"""
    dG = sum(nu * DB.G(sp, T) for sp, nu in REACTION.items())
    return np.exp(-dG / (R * T))


def gas_state(T, hcl_o2_ratio=2.0, P=1.0, extent=None) -> GasState:
    """供給 HCl/O2 比に対する気相状態。

    Parameters
    ----------
    extent : None なら平衡まで進んだ状態。数値ならその反応進行度。
    """
    n_hcl0, n_o20 = float(hcl_o2_ratio), 1.0
    xi_max = min(n_hcl0 / 4, n_o20)

    def partials(xi):
        n = {
            "HCl(g)": n_hcl0 - 4 * xi,
            "O2(g)": n_o20 - xi,
            "Cl2(g)": 2 * xi,
            "H2O(g)": 2 * xi,
        }
        tot = sum(n.values())
        return {k: max(v / tot * P, 1e-30) for k, v in n.items()}

    if extent is None:
        K = equilibrium_constant(T)

        def residual(xi):
            p = partials(xi)
            lnQ = (
                2 * np.log(p["Cl2(g)"]) + 2 * np.log(p["H2O(g)"])
                - 4 * np.log(p["HCl(g)"]) - np.log(p["O2(g)"])
            )
            return lnQ - np.log(K)

        xi = brentq(residual, 1e-12, xi_max * (1 - 1e-12))
    else:
        xi = float(extent)

    p = partials(xi)
    return GasState(
        T=T, p_HCl=p["HCl(g)"], p_O2=p["O2(g)"], p_Cl2=p["Cl2(g)"],
        p_H2O=p["H2O(g)"], conversion=4 * xi / n_hcl0, extent=xi,
    )


def operating_line(T, hcl_o2_ratio=2.0, P=1.0, n=200) -> list[GasState]:
    """入口（転化率ゼロ）から平衡までの軌跡。安定領域図に重ねて使う。"""
    xi_eq = gas_state(T, hcl_o2_ratio, P).extent
    return [
        gas_state(T, hcl_o2_ratio, P, extent=x)
        for x in np.linspace(1e-4 * xi_eq, xi_eq * (1 - 1e-9), n)
    ]
