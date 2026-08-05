"""融液上の Cu 蒸気圧と触媒寿命。

塩化物系 Deacon 触媒の主たる失活経路は Cu の揮発であり、
KCl / LnCl3 を加える主目的はこれを抑えることにある。

    p(Cu 種) = K(T) * a(CuCl or CuCl2)^n

なので必要なのは気相種の生成自由エネルギー（data.py）と
融液中の活量（melt.py）の二つ。
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .data import DB
from .melt import Melt
from .species import R

#: 気相種あたりの Cu 原子数
CU_VAPOUR_SPECIES = {"CuCl(g)": ("CuCl", 1), "Cu3Cl3(g)": ("CuCl", 3),
                     "CuCl2(g)": ("CuCl2", 1)}

R_ATM = 0.0820574  # L·atm/(mol·K)
M_CU = 63.546  # g/mol


def partial_pressures(melt: Melt, T) -> dict[str, float]:
    """融液上の Cu 含有気相種の分圧 [atm]。"""
    out = {}
    for gas, (salt, n) in CU_VAPOUR_SPECIES.items():
        liquid = "CuCl(s)" if salt == "CuCl" else "CuCl2(s)"
        dG = DB.G(gas, T) - n * DB[liquid].G_supercooled_liquid(T)
        out[gas] = melt.activity(salt, T) ** n * np.exp(-dG / (R * T))
    return out


def cu_vapour_fraction(melt: Melt, T, P=1.0) -> float:
    """気相 1 mol あたりに持ち去られる Cu のモル数。"""
    p = partial_pressures(melt, T)
    return sum(p[g] * n for g, (_, n) in CU_VAPOUR_SPECIES.items()) / P


@dataclass(frozen=True)
class ReactorSpec:
    """揮発速度を寿命に換算するための反応器条件。"""

    ghsv: float = 6000.0  # L/(kg-cat·h)
    cu_loading_wt_pct: float = 10.0
    P: float = 1.0

    def gas_flow(self, T) -> float:
        """mol/(kg-cat·h)"""
        return self.ghsv * self.P / (R_ATM * T)

    @property
    def cu_inventory(self) -> float:
        """mol Cu/kg-cat"""
        return self.cu_loading_wt_pct * 10.0 / M_CU


def lifetime(melt: Melt, T, spec: ReactorSpec | None = None, frac_lost=0.5):
    """Cu を frac_lost だけ失うまでの時間 [h] と損失速度 [mg/(kg-cat·h)]。"""
    spec = spec or ReactorSpec()
    rate = cu_vapour_fraction(melt, T, spec.P) * spec.gas_flow(T)  # mol/(kg·h)
    if rate <= 0:
        return np.inf, 0.0
    return frac_lost * spec.cu_inventory / rate, rate * M_CU * 1000


def required_activity_coefficient(observed_lifetime_h, melt: Melt, T,
                                  spec: ReactorSpec | None = None, frac_lost=0.5):
    """逆問題: 実測寿命を再現するのに必要な活量係数を逆算する。

    理想 Temkin が実測を何桁外すかが、そのまま
    「クロロ銅酸錯体による安定化 + Cu3Cl3(g) データ誤差」の大きさになる。
    両者を分離するには文献の平衡塩素圧データが要る。
    """
    spec = spec or ReactorSpec()
    rate_needed = frac_lost * spec.cu_inventory / observed_lifetime_h
    flux_needed = rate_needed / spec.gas_flow(T)
    flux_ideal = cu_vapour_fraction(melt, T, spec.P)
    ratio = flux_needed / flux_ideal
    return {
        "flux_ideal": flux_ideal,
        "flux_needed": flux_needed,
        "ratio": ratio,
        "orders_of_magnitude": float(np.log10(ratio)),
        "gamma_if_trimer_dominant": ratio ** (1 / 3),
        "gamma_if_monomer_dominant": ratio,
    }
