# deacon-thermo

Cu-K-Ln/γ-Al₂O₃ 系 HCl 酸化（Deacon）触媒の熱力学解析。

反応条件下での凝縮相の安定性、Cu の揮発による失活速度、
ランタノイド置換による変化を計算する。

依存は numpy と scipy のみ。FactSage も pycalphad も不要。

> **先に [CLAUDE.md](CLAUDE.md) を読むこと。** この系は塩化物三元系ではなく
> Cu-K-Ln-Cl-O-H 系であるという前提を外すと、計算は動くが意味のない結果が出る。

## インストール

```bash
pip install -e ".[plot,analysis,dev]"
pytest
```

## 使い方

```python
from deacon_thermo import gas_state, hydrolysis_ranking, redox_split, lifetime

# 380 C, HCl/O2 = 2:1 の平衡気相
gas = gas_state(T=653.15, hcl_o2_ratio=2.0)
print(gas.conversion, gas.p_Cl2)

# Ln を塩化物として残りやすい順に並べる
for d in hydrolysis_ranking(gas):
    print(f"{d.element:3s} margin={d.chloride_margin:+7.1f} kJ/mol  r={d.ionic_radius}")

# 融液上の Cu 揮発と寿命
_, melt = redox_split(cu_total=1.0, p_Cl2=gas.p_Cl2, T=653.15,
                      diluents={"KCl": 1.0, "SmCl3": 0.3})
t_half, rate = lifetime(melt, 653.15)
```

例は `examples/` を参照。

```bash
python examples/01_stability_diagram.py
python examples/02_lanthanide_survey.py
```

## データの信頼度

```python
from deacon_thermo import DB
print(DB.report())
```

`confidence == ESTIMATE` の種は結論に効くなら一次資料で置換すること。
現時点で SmOCl（および他の LnOCl）と Cu₃Cl₃(g) が最重要。

## 既知の制約

- 融液の活量モデルは理想 Temkin が既定で、錯体形成による安定化を含まない。
  Cu 蒸気圧を約 2 桁過大評価する（`test_ideal_model_reproduces_observed_lifetime`
  が xfail として残してある）。`fit_interactions()` に平衡塩素圧データを
  入れて較正すること。
- 液相線は未実装。必要になったら `melt.ActivityModel` のサブクラスとして
  pycalphad の MQMQA を接続する。
- γ-Al₂O₃ 担体との反応・K/Cu の偏析は熱力学計算の外にある。
