"""Cu-K-Ln/Al2O3 系 Deacon 触媒の熱力学解析。

この系は塩化物三元系ではなく Cu-K-Ln-Cl-O-H 系であり、
CuCl2 = CuCl + 1/2 Cl2 は障害ではなく反応機構そのものである。
背景と設計判断は CLAUDE.md を参照。

典型的な使い方::

    from deacon_thermo import gas_state, hydrolysis_ranking

    gas = gas_state(T=653.15, hcl_o2_ratio=2.0)
    for d in hydrolysis_ranking(gas):
        print(d.element, d.chloride_margin)
"""

from .assemblage import (
    Assemblage,
    CationGrid,
    assemblage_at,
    cation_grid,
    default_candidates,
    stable_assemblage,
)
from .data import CALIBRATED_INTERACTIONS, DB, LANTHANIDES
from .gas import GasState, equilibrium_constant, gas_state, operating_line
from .lanthanides import (
    REFERENCE_CATIONS,
    LanthanideDescriptors,
    chloride_fraction,
    descriptors,
    hydrolysis_ranking,
    k3lncl6_stability_limit,
    oxychloride_params,
    radius_controls,
    survey,
)
from .melt import (
    ClObservation,
    IdealTemkin,
    Melt,
    RegularSolution,
    calibrated_model,
    equilibrium_p_cl2,
    fit_interactions,
    redox_K,
    redox_split,
)
from .sensitivity import (
    MarginSweep,
    chloride_margin_at,
    flip_threshold,
    flip_thresholds,
    sweep_margins,
)
from .species import Confidence, Species
from .stability import (
    copper_phases,
    dHf_oxychloride_threshold,
    hydrolysis_dG,
    hydrolysis_K,
    hydrolysis_margin,
    lanthanide_phases,
    stability_map,
    stable_chloride,
)
from .volatility import (
    ReactorSpec,
    cu_vapour_fraction,
    lifetime,
    partial_pressures,
    required_activity_coefficient,
)

__version__ = "0.1.0"

__all__ = [
    "DB", "LANTHANIDES", "CALIBRATED_INTERACTIONS", "Species", "Confidence",
    "GasState", "gas_state", "operating_line", "equilibrium_constant",
    "lanthanide_phases", "copper_phases", "stability_map",
    "hydrolysis_dG", "hydrolysis_K", "stable_chloride", "hydrolysis_margin",
    "dHf_oxychloride_threshold",
    "Melt", "IdealTemkin", "RegularSolution", "calibrated_model",
    "redox_split", "fit_interactions",
    "ClObservation", "equilibrium_p_cl2", "redox_K",
    "partial_pressures", "cu_vapour_fraction", "lifetime", "ReactorSpec",
    "required_activity_coefficient",
    "LanthanideDescriptors", "descriptors", "survey", "hydrolysis_ranking",
    "radius_controls", "REFERENCE_CATIONS", "oxychloride_params",
    "k3lncl6_stability_limit", "chloride_fraction",
    "MarginSweep", "chloride_margin_at", "sweep_margins",
    "flip_threshold", "flip_thresholds",
    "Assemblage", "CationGrid", "assemblage_at", "stable_assemblage",
    "cation_grid", "default_candidates",
]
