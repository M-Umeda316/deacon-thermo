"""安定領域図と操業線。

    python examples/01_stability_diagram.py
"""

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from deacon_thermo import (  # noqa: E402
    DB,
    copper_phases,
    dHf_oxychloride_threshold,
    gas_state,
    hydrolysis_margin,
    lanthanide_phases,
    operating_line,
    stability_map,
)

T = 653.15
LN = "Sm"


def main():
    gas = gas_state(T, hcl_o2_ratio=2.0)
    line = operating_line(T, hcl_o2_ratio=2.0)

    print(f"T = {T - 273.15:.0f} C,  HCl/O2 = 2:1,  P = 1 atm")
    print(f"  平衡 HCl 転化率 = {gas.conversion * 100:.1f} %")
    for k, v in gas.as_dict().items():
        print(f"    p({k:7s}) = {v:.4f} atm")

    print(f"\n{LN}Cl3 + H2O = {LN}OCl + 2HCl")
    for label, g in [("入口", line[0]), ("出口(平衡)", gas)]:
        m = hydrolysis_margin(LN, g)
        thr = dHf_oxychloride_threshold(LN, g)
        verdict = f"{LN}Cl3" if m > 0 else f"{LN}OCl"
        print(f"  {label:10s} 余裕 {m:+8.1f} kJ/mol -> {verdict:8s}"
              f"  (dHf({LN}OCl) < {thr:.1f} で反転)")
    print(f"  現在の仮値: {DB[f'{LN}OCl(s)'].dHf298:.1f} kJ/mol")
    print("\n" + DB.report())

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
    for ax, (phases, title) in zip(
        axes, [(lanthanide_phases(LN), f"{LN}-Cl-O"), (copper_phases(), "Cu-Cl-O")], strict=True
    ):
        smap = stability_map(phases, T)
        cmap = plt.get_cmap("Pastel1")
        ax.contourf(smap.log_pO2, smap.log_pCl2, smap.index,
                    levels=np.arange(-0.5, len(smap.names)),
                    colors=[cmap(i) for i in range(len(smap.names))])
        ax.contour(smap.log_pO2, smap.log_pCl2, smap.index,
                   levels=np.arange(0.5, len(smap.names)), colors="k", linewidths=0.9)
        for i, nm in enumerate(smap.names):
            mask = smap.index == i
            if mask.sum() > 40:
                yi, xi = np.where(mask)
                ax.text(smap.log_pO2[int(np.median(xi))],
                        smap.log_pCl2[int(np.median(yi))],
                        nm, ha="center", va="center", fontsize=11, weight="bold")

        ax.plot([np.log10(g.p_O2) for g in line], [np.log10(g.p_Cl2) for g in line],
                "r-", lw=2.4, zorder=5, label="operating line")
        ax.plot(np.log10(gas.p_O2), np.log10(gas.p_Cl2), "r*", ms=17, zorder=6,
                label="equilibrium")
        ax.set_xlabel("log $p_{\\mathrm{O_2}}$ / atm")
        ax.set_ylabel("log $p_{\\mathrm{Cl_2}}$ / atm")
        ax.set_title(f"{title}  ({T - 273.15:.0f} $^\\circ$C)")
        ax.legend(loc="lower right", fontsize=9)

    fig.tight_layout()
    fig.savefig("stability_380C.png", dpi=160)
    print("\n-> stability_380C.png")


if __name__ == "__main__":
    main()
