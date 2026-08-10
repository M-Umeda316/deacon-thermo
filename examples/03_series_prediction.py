"""Ln 系列の統一予測: 相・複塩・Cu 揮発寿命を 1 枚の表にする。

    python examples/03_series_prediction.py

3 つの問い（CLAUDE.md）を同じ条件で並べる:
  1. Ln は塩化物のままか      -> chloride_margin と安定相集合
  2. Cu はどれだけ揮発するか  -> 較正済み融液モデルでの半減時間
  3. Ln を変えるとどう動くか  -> 系列を縦に並べて比較

寿命だけ 613 K / p(Cl2) = 0.32 atm で評価するのは、参照系（Cu-K-La、9600 h 安定）
の報告条件に合わせて実測と突き合わせられるようにするため。相の判定は
反応条件（380 C）の平衡出口気相で行う。
"""

from deacon_thermo import (
    LANTHANIDES,
    ReactorSpec,
    gas_state,
    hydrolysis_ranking,
    k3lncl6_stability_limit,
    lifetime,
    redox_split,
    stable_assemblage,
)
from deacon_thermo.data import DB
from deacon_thermo.melt import RegularSolution, calibrated_model

T_GAS = 653.15  # 380 C、相の判定
T_VOL = 613.15  # 340 C、参照系の寿命報告条件
P_CL2_VOL = 0.32  # 同上
CATIONS = {"Cu": 0.4, "K": 0.35, "Ln": 0.25}
DILUENTS = {"KCl": 1.0, "LnCl3": 0.3}

#: W(Cu 塩化物, LnCl3) は未較正。感度帯として振る幅 [J/mol]
W_CU_LN_BAND = (-10000.0, 0.0, +10000.0)


def melt_model(ln: str, w_cu_ln: float) -> RegularSolution:
    """較正済みモデルに W(CuCl2, LnCl3) を足したもの。"""
    model = calibrated_model()
    if w_cu_ln:
        model.interactions[("CuCl2", f"{ln}Cl3")] = w_cu_ln
    return model


def half_lives(ln: str, spec: ReactorSpec | None = None) -> list[float]:
    """W(Cu-Ln) を感度帯で振ったときの Cu 半減時間 [h]。"""
    diluents = {"KCl": DILUENTS["KCl"], f"{ln}Cl3": DILUENTS["LnCl3"]}
    out = []
    for w in W_CU_LN_BAND:
        _, melt = redox_split(1.0, P_CL2_VOL, T_VOL, diluents, melt_model(ln, w))
        out.append(lifetime(melt, T_VOL, spec)[0])
    return out


def assemblage_summary(ln: str, gas) -> str:
    """相集合を「相名 (Ln 分配%)」の短い文字列にする。"""
    asm = stable_assemblage({"Cu": CATIONS["Cu"], "K": CATIONS["K"], ln: CATIONS["Ln"]}, gas)
    parts = []
    for name, amount in sorted(asm.phases.items(), key=lambda kv: -kv[1]):
        nu = DB[name].elements.get(ln, 0.0)
        short = name.removesuffix("(s)")
        if nu:
            parts.append(f"{short} {100 * amount * nu / CATIONS['Ln']:.0f}%")
        else:
            parts.append(short)
    return " + ".join(parts)


def main():
    gas = gas_state(T_GAS, hcl_o2_ratio=2.0)
    spec = ReactorSpec()
    margins = {d.element: d.chloride_margin for d in hydrolysis_ranking(gas)}

    print(f"Ln 系列の統一予測  (相: {T_GAS - 273.15:.0f} C, HCl/O2 = 2:1, 平衡出口 "
          f"p(Cl2) = {gas.p_Cl2:.3f} atm)")
    print(f"  カチオン組成 Cu:K:Ln = "
          f"{CATIONS['Cu']:.2f}:{CATIONS['K']:.2f}:{CATIONS['Ln']:.2f}")
    print(f"  寿命: {T_VOL - 273.15:.0f} C, p(Cl2) = {P_CL2_VOL} atm, "
          f"融液 CuCl+CuCl2 1 : KCl {DILUENTS['KCl']} : LnCl3 {DILUENTS['LnCl3']}, "
          f"Cu 半減 (GHSV {spec.ghsv:.0f})\n")

    header = (f"{'Ln':4s} {'余裕/kJ':>8s}  {'安定相集合 (Ln 分配)':<34s} "
              f"{'K3LnCl6':>9s}  {'t_half/h  W(Cu-Ln)= -10k / 0 / +10k':>36s}")
    print(header)
    print("-" * 105)

    for ln in sorted(LANTHANIDES, key=lambda e: -margins[e]):
        bound = k3lncl6_stability_limit(ln)
        bound_txt = f"{bound[0]:.0f}{bound[1]}" if bound else "-"
        lo, mid, hi = half_lives(ln, spec)
        print(f"{ln:4s} {margins[ln]:+8.1f}  {assemblage_summary(ln, gas):<34s} "
              f"{bound_txt:>9s}  {lo:11.3g}{mid:12.3g}{hi:12.3g}")

    print("\nK3LnCl6 列: 末尾 L = L/H 転移温度（低温形が別にある）、"
          "S = 合成 dG=0（それ以下では二元塩化物に分解）")
    print("\n注意（この表を読むときの前提）")
    print("  - W(Cu-LnCl3) は未較正。左右の列は ±10 kJ/mol を仮に振った帯であって")
    print("    誤差棒ではない。Ruthven & Kenney の四元系データが入るまで暫定。")
    print("  - 較正が効いているのは Cu-K 系（KCl 30 mol%）だけなので、LnCl3 は今のところ")
    print("    理想希釈剤にすぎず、寿命 3 列は Ln に依らず同じ値になる。ここに系列差を")
    print("    出すには W(Cu-LnCl3) の較正が要る、というのがこの列の読み方。")
    print("  - 反応条件の支配蒸気種 CuCl2(g) は data.py で ESTIMATE。寿命の絶対値は")
    print("    このデータに直結するので、桁の議論にしか使わないこと。")
    print("  - assemblage の候補相は固体のみ。380 C では Cu-K 塩化物は実際には融体で、")
    print("    K2CuCl3 はその代理でしかない（液相線は未実装）。")


if __name__ == "__main__":
    main()
