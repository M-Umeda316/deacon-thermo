"""熱力学量の基本クラス。

Maier-Kelley 型の熱容量から H(T), S(T), G(T) を組み立てる。

    Cp(T) = a + b*1e-3*T + c*1e5/T^2   [J/mol/K]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

R = 8.314462618  # J/mol/K
T_REF = 298.15  # K


class Confidence(str, Enum):
    """データの信頼度。ESTIMATE のものは結論に効くなら必ず一次資料で置換する。"""

    GOOD = "good"  # NIST-JANAF / Barin 由来、そのまま使える
    FAIR = "fair"  # 文献値だが版・出典により差がある
    ESTIMATE = "estimate"  # 推定値。要検証


@dataclass(frozen=True)
class Species:
    """単一の化学種。

    Attributes
    ----------
    dHf298 : 標準生成エンタルピー [kJ/mol]
    S298   : 標準エントロピー [J/mol/K]
    cp     : Maier-Kelley 係数 (a, b, c)
    elements : 元素組成 {'Sm': 1, 'Cl': 3}
    transitions : [(転移温度 K, 転移エンタルピー kJ/mol), ...]
    """

    name: str
    phase: str  # 's' | 'l' | 'g'
    dHf298: float
    S298: float
    cp: tuple[float, float, float]
    elements: dict[str, float]
    confidence: Confidence = Confidence.FAIR
    source: str = ""
    transitions: tuple[tuple[float, float], ...] = field(default_factory=tuple)

    # -- 熱力学量 ----------------------------------------------------------
    def _cp_int(self, T):
        a, b, c = self.cp
        b, c = b * 1e-3, c * 1e5
        return a * (T - T_REF) + b / 2 * (T**2 - T_REF**2) - c * (1 / T - 1 / T_REF)

    def _cp_over_T_int(self, T):
        a, b, c = self.cp
        b, c = b * 1e-3, c * 1e5
        return (
            a * np.log(T / T_REF)
            + b * (T - T_REF)
            - c / 2 * (1 / T**2 - 1 / T_REF**2)
        )

    def H(self, T):
        """H(T) [J/mol]、298 K の元素基準。"""
        h = self.dHf298 * 1000 + self._cp_int(T)
        for T_tr, dH_tr in self.transitions:
            h = h + np.where(T >= T_tr, dH_tr * 1000, 0.0)
        return h

    def S(self, T):
        """S(T) [J/mol/K]。"""
        s = self.S298 + self._cp_over_T_int(T)
        for T_tr, dH_tr in self.transitions:
            s = s + np.where(T >= T_tr, dH_tr * 1000 / T_tr, 0.0)
        return s

    def G(self, T):
        """G(T) [J/mol]。"""
        return self.H(T) - T * self.S(T)

    def G_supercooled_liquid(self, T):
        """過冷却液体を基準にした G。融液成分の標準状態に使う。

        G_l(T) = G_s(T) + dH_fus * (1 - T/T_fus)
        transitions の最初の要素を融解とみなす。
        """
        if not self.transitions:
            raise ValueError(f"{self.name} に融解データがない")
        T_fus, dH_fus = self.transitions[0]
        # 固体としての G（転移の段差を含めない）
        G_s = (self.dHf298 * 1000 + self._cp_int(T)) - T * (
            self.S298 + self._cp_over_T_int(T)
        )
        return np.where(T < T_fus, G_s + dH_fus * 1000 * (1 - T / T_fus), G_s)


class Database(dict):
    """化学種のレジストリ。"""

    def add(self, species: Species) -> Species:
        if species.name in self:
            raise KeyError(f"{species.name} は既に登録済み")
        self[species.name] = species
        return species

    def G(self, name: str, T):
        return self[name].G(T)

    def needs_verification(self) -> list[Species]:
        """信頼度 ESTIMATE の種を列挙する。"""
        return [s for s in self.values() if s.confidence is Confidence.ESTIMATE]

    def report(self) -> str:
        lines = ["要検証のデータ:"]
        for s in self.needs_verification():
            lines.append(f"  {s.name:14s} dHf298={s.dHf298:9.1f} kJ/mol  {s.source}")
        return "\n".join(lines)
