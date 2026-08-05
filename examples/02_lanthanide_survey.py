"""ランタノイド系列の記述子表と実験計画への落とし込み。

    python examples/02_lanthanide_survey.py
"""

from deacon_thermo import (
    ReactorSpec,
    cu_vapour_fraction,
    gas_state,
    hydrolysis_ranking,
    lifetime,
    radius_controls,
    redox_split,
)

T = 653.15


def main():
    gas = gas_state(T, hcl_o2_ratio=2.0)

    print(f"Ln 系列の比較  ({T - 273.15:.0f} C, HCl/O2 = 2:1, 平衡出口)")
    print("塩化物として残りやすい順。系列内の差は絶対値より信頼できる。\n")
    print(f"{'Ln':4s} {'r/Å':>7s} {'酸化還元':>8s} {'余裕/kJ':>9s} "
          f"{'閾値dHf':>10s} {'LnCl3 mp/K':>11s}")
    print("-" * 56)
    for d in hydrolysis_ranking(gas):
        flag = "+".join(f"{s}" for s in d.accessible_states) if d.accessible_states else "-"
        print(f"{d.element:4s} {d.ionic_radius:7.3f} {flag:>8s} "
              f"{d.chloride_margin:+9.1f} {d.dHf_threshold:10.1f} {d.melting_point:11.0f}")

    print("\n[軸A] 酸化還元活性な Ln は機構が変わる。")
    print("      Ce, Pr は酸化雰囲気でオキシ塩化物にならず CeO2 / PrO2 に行く")
    print("      （このモジュールの候補相には未実装 — 別扱いが要る）")

    print("\n[半径対照ペア] 4f の寄与を切り分けられる組み合わせ:")
    for a, b in radius_controls():
        print(f"      {a} / {b}")
    print("      特に Y / Ho: Y は 4f を持たないので決定実験になる")

    print("\n[Cu 揮発への希釈効果] Ln の種類ではなく量で効く部分")
    spec = ReactorSpec()
    for label, dil in [
        ("希釈なし", {}),
        ("KCl 1", {"KCl": 1.0}),
        ("KCl 1 + LnCl3 0.3", {"KCl": 1.0, "LnCl3": 0.3}),
        ("KCl 2 + LnCl3 1.0", {"KCl": 2.0, "LnCl3": 1.0}),
    ]:
        _, melt = redox_split(1.0, gas.p_Cl2, T, dil)
        t_half, rate = lifetime(melt, T, spec)
        print(f"  {label:20s} Cu 蒸気分率 {cu_vapour_fraction(melt, T):.3e}  "
              f"損失 {rate:.2e} mg/kg/h  半減 {t_half:.2e} h")
    print("\n  注: 理想 Temkin なので絶対値は約2桁過大。比較のみに使うこと。")


if __name__ == "__main__":
    main()
