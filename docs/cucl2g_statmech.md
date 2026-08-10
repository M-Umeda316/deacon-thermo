# CuCl₂(g) の熱力学関数を統計力学で拘束する

作成日: 2026-08-10 / 計算環境: conda `deacon-thermo`
対象: `data.py` の `CuCl2(g)`（現行 ΔHf₂₉₈ = -43.0 kJ/mol, S°₂₉₈ = 278.0 J/mol/K,
Cp = MK(57.0, 0.5, -0.5), confidence = ESTIMATE）

Cu₃Cl₃(g) を JANAF 値に差し替えた結果、反応条件（613-653 K、p(Cl₂) ~ 0.3 atm）で
Cu を運ぶ気相種は **CuCl₂(g) が 100%** になった（本文 §4 の分圧内訳を参照）。
したがって揮発の結論は CuCl₂(g) の ΔG(T) にほぼ完全にぶら下がっている。
本ノートは、そのうち **S°₂₉₈ と Cp(T) を分光定数から確定**し、
**ΔHf₂₉₈ については実測寿命から下限を与える**ことを目的とする。

結論を先に書く:

| 量 | 現行暫定 | 本ノート | 判定 |
|---|---|---|---|
| S°₂₉₈ [J/mol/K] | 278.0 (EST) | **285.5 ± 4** (統計力学) | 差 +7.5 → 置換すべき |
| Cp₂₉₈ [J/mol/K] | 56.6 | **60.6** | 差 +4.0 |
| Cp₆₅₃ [J/mol/K] | 57.2 | **64.2** | 差 +7.0 |
| p(CuCl₂) @653 K | — | **×2.87**（同一 ΔHf で） | 揮発は現行見積りより速い |
| ΔHf₂₉₈ [kJ/mol] | -43.0 (EST) | **≥ -52.0 〜 -49.3**（下限、実験条件 GHSV 450） | **現行値 -43.0 は整合**（検証追記参照） |

---

## 検証追記(2026-08-10、親エージェントによる訂正)

本ノート初版の ΔHf 下限(≥ -38.8 〜 -36.1)と「現行値 -43.0 は不適格」の判定は、
逆算に **ReactorSpec 既定の GHSV 6000 L/kg/h** を使っていた。しかし基準にした
実測「9600 h 安定」(Feng ら 2015, Appl. Catal. B 164, 483) の実条件は
**GHSV 450 L/kg/h**(当プロジェクトの文献調査で記録済み)であり、
下限は RT ln(6000/450) = 13.2 kJ/mol 緩む:

- 実験条件 (GHSV 450) での下限: **ΔHf ≥ -52.0 (理想 Temkin) / -49.3 (較正モデル)**
- 現行値 -43.0 での予測 t_half (新 S/Cp、GHSV 450): 5.7×10⁴ h (理想) /
  3.3×10⁴ h (較正) — **実測 ≥9600 h と矛盾しない**

したがって **-43.0 は棄却されず EST のまま保持**し、S°/Cp のみ本ノートの
統計力学値へ差し替えた(data.py 反映済み)。§4 の初版の数値は GHSV 6000 の
設計点としての計算として読むこと。S°₂₉₈ = 285.46 と各寄与は親エージェントが
独立再計算して 0.01 J/mol/K まで一致を確認済み。

---

## 1. 採用した分子定数

CuCl₂ は遊離分子では直線 D∞h、電子基底状態は X ²Π_g（スピン軌道分裂は逆転、
²Π_3/2g が下）。すべて実測値で、推定で埋めた定数は無い。

| 量 | 採用値 | 媒体/方法 | 出典 | URL | 取得日 | 信頼度 |
|---|---|---|---|---|---|---|
| r₀(Cu-Cl) | 203.52 pm | 高分解能 LIF + 同位体置換 | Beattie, Brown, Crozet, Ross, Yiannopoulou, *Inorg. Chem.* **36** (1997) 3207-3208, "On the Geometry of the CuCl₂ Molecule" | https://pubs.acs.org/doi/10.1021/ic961509r （本文は有料。値は NIST WebBook / 公開書誌の記載から） | 2026-08-10 | GOOD |
| ν₁ (σ_g⁺, 対称伸縮) | 369.36 ± 0.04 cm⁻¹ | 気相, LF | NIST WebBook, CuCl₂ 振動/電子準位表（Barnes & Brown 2007。同表に ω₁ 相当の 371.69 も併記） | https://webbook.nist.gov/cgi/cbook.cgi?ID=B1001128&Mask=800 | 2026-08-10 | GOOD |
| ν₂ (π_u, 変角, 2 重縮退) | 95.81 cm⁻¹ | 気相, LF（Ne 97.0 / Ar 98.6） | 同上（Crozet, Ross et al. 1995） | 同上 | 2026-08-10 | **FAIR**（§5 参照） |
| ν₃ (σ_u⁺, 反対称伸縮) | 522.18 cm⁻¹ | 気相（同表に 525.90 も。Ne 522.3 / Ar 514.4, IR） | 同上（Barnes & Brown 2007 / マトリックス IR） | 同上 | 2026-08-10 | GOOD |
| 対称数 σ | 2 | D∞h | — | — | — | — |
| 分子量 | 134.452 g/mol | 標準原子量 | NIST WebBook | 同上 | 2026-08-10 | GOOD |

### 電子状態（分配関数に入れたもの）

Ω 成分ごとに縮重 2（±Ω）として数える。

| 状態 | E [cm⁻¹] | g | 出典 |
|---|---|---|---|
| X ²Π_3/2 g | 0 | 2 | 基底 |
| ²Π_1/2 g | **482.9** | 2 | Hodges & Ross, *J. Chem. Phys.* **127** (2007) 024309（気相 LIF、A = -482.9 cm⁻¹）。NIST WebBook は Ne マトリックス値 474.1、Ar 303.6 を併記 |
| ²Σ_g⁺ | **1910.9** (Ne) | 2 | Lorenz, Caspary et al. 1997、NIST WebBook 記載。気相値は未測定（Ar では 1616.5） |
| ²Δ_5/2 g | 6877 (気相) | 2 | Wang et al. 2001 ほか、NIST WebBook |
| ²Δ_3/2 g | 8700 ± 900 (気相) | 2 | 同上 |
| ²Σ_u⁺ / ²Π_u ほか | ≥ 15500 | — | 653 K で寄与 < 10⁻¹⁰、無視 |

出典 URL（電子状態）:
- NIST Chemistry WebBook, CuCl₂ 振動・電子エネルギー準位: https://webbook.nist.gov/cgi/cbook.cgi?ID=B1001128&Mask=800 （2026-08-10 取得）
- Hodges & Ross (2007) 書誌: https://pubs.aip.org/aip/jcp/article/127/2/024309/919449 （本文 403、要旨から A = -482.9 cm⁻¹ を確認、2026-08-10）

**注意（未決着）**: ²Σ_g⁺ の位置は理論と実験で大きく食い違う。CASSCF+ACPF 系の
ベンチマーク計算では X ²Π_g → ²Σ_g⁺ の垂直遷移が **99 cm⁻¹** と出ており、
文献全体では -1856 〜 +5887 cm⁻¹ の幅がある
（Ab initio study on the spectroscopy of CuCl₂ I/II、https://pubmed.ncbi.nlm.nih.gov/15268475/ 、
https://pubmed.ncbi.nlm.nih.gov/15945683/ 、2026-08-10 取得。OA 全文 https://hal.science/hal-00003672
は取得時 Anubis で遮断されアクセスできず）。
本ノートは実測（Ne マトリックス 1910.9）を採ったが、感度は §5 に出した。

---

## 2. 統計力学計算

剛体回転子（直線, σ = 2）+ 調和振動子 + 電子分配関数。標準状態は 1 atm
（JANAF 準拠）。原子量は標準原子量（同位体平均）を使う — 熱力学関数の慣行に合わせる。

慣性モーメント: I = 2 m_Cl r² = 48.769 × 10⁻⁴⁶ kg m²（B = 0.05740 cm⁻¹）

### 2.1 コードの検証（CuCl(g)）

同じコードで CuCl(g)（r = 2.0512 Å, ω = 415.3 cm⁻¹, 基底 ¹Σ⁺, σ = 1）を計算:

| 量 | 本計算 | 参照 (JANAF / CCCBDB) | 差 |
|---|---|---|---|
| S°₂₉₈ | **236.96** | 237.20 | **-0.24 J/mol/K** |
| Cp₂₉₈ | 35.11 | 35.26 | -0.15 J/mol/K |

内訳: 並進 166.05 / 回転 67.10 / 振動 3.80 / 電子 0.00。
要求の ±1 J/mol/K を満たす。**検証成功**。
（回転定数は同位体平均質量で 0.17606 cm⁻¹ となり、CCCBDB の ⁶³Cu³⁵Cl 値 0.17825 とは
同位体の取り方の違いだけずれる。熱力学関数側では平均質量が正しい。）

参照値出典: NIST CCCBDB, CuCl 実験データ
https://cccbdb.nist.gov/exp2x.asp?casno=7758896 （2026-08-10 取得）

### 2.2 CuCl₂(g) の S°₂₉₈

| 寄与 | S [J/mol/K] |
|---|---|
| 並進 | 169.87 |
| 回転 | 70.66 |
| 振動 | 36.66 |
| 電子 | 8.26 |
| **合計 S°₂₉₈** | **285.46** |

振動の 36.66 のうち **29.6 は ν₂（95.81 cm⁻¹, 2 重縮退）** から来る。
電子の 8.26 は主に基底の 2 重縮重（R ln2 = 5.76）と ²Π_1/2 の熱励起。

現行暫定値 278.0 との差: **+7.46 J/mol/K**。

（調和振動数版 ω = (371.69, 95.81, 525.90) を使っても 285.38 で、差は 0.08。
伸縮の調和/基本振動数の区別はここでは効かない。）

### 2.3 Cp(T)

| T [K] | Cp [J/mol/K] | MK フィット | 現行暫定 | 差 |
|---|---|---|---|---|
| 298.15 | 60.592 | 60.58 | 56.59 | +4.00 |
| 300 | 60.646 | | | |
| 400 | 62.516 | | | |
| 500 | 63.404 | | | |
| 600 | 63.950 | | | |
| **653.15** | **64.161** | | 57.21 | **+6.95** |
| 700 | 64.311 | | | |
| 800 | 64.536 | | | |
| 900 | 64.659 | | | |
| 1000 | 64.710 | 64.68 | 57.45 | +7.26 |

Maier-Kelley フィット（Cp = a + b·10⁻³·T + c·10⁵/T²、298-1000 K）:

```
(a, b, c) = (64.845, 0.3103, -3.8724)     最大残差 0.058 J/mol/K
```

参考: 300-1000 K のみだと (64.858, 0.2973, -3.8905)、298-1200 K だと
(65.133, -0.0533, -4.0452)（最大残差 0.087）。反応温度域では実質同じ。

古典極限 Cp = 4R + 4R = ... ではなく、直線 3 原子分子の上限は
7/2 R + 4R = 62.4 J/mol/K。計算値が 64 を超えるのは電子項（Schottky）の寄与。

---

## 3. p(CuCl₂) への影響

ΔHf を -43.0 に固定したまま S° と Cp を上の計算値に置き換えると、
CuCl₂(l) → CuCl₂(g) の ΔG が下がり、蒸気圧が上がる:

| T [K] | ΔG(g) の変化 | p(CuCl₂) 比 |
|---|---|---|
| 613.15 | -5254.7 J/mol | **×2.80** |
| 653.15 | -5727.8 J/mol | **×2.87** |

つまり **現行の DB は Cu の揮発をおよそ 3 倍過小評価している**（ΔHf が正しい場合）。
寿命はその逆数で効くので、参照触媒に対する予測寿命は 1.19×10⁴ h → 4.25×10³ h
（IdealTemkin）に短くなる。

---

## 4. ΔHf₂₉₈ の下限（実測寿命からの逆算）

### 条件

- Cu-K-La/γ-Al₂O₃ 参照触媒: **340 ℃ = 613.15 K**、HCl 転化率 78% を **9600 h 維持**
- `redox_split(1.0, 0.32, 613.15, {'KCl': 1.0, 'LaCl3': 0.3}, model)` の融液
- `ReactorSpec()` 既定（GHSV 6000 L/(kg·h)、Cu 10 wt%、1 atm）
- 「9600 h 安定」を t_half の**下限**とみなす → lifetime ≥ 9600 h を要求
- CuCl(g) / Cu₃Cl₃(g) は現行 DB のまま

### 融液の状態と分圧内訳（ΔHf = -43.0、新 S/Cp）

| モデル | f(Cu-II) | a(CuCl) | a(CuCl₂) | p(CuCl) | p(Cu₃Cl₃) | p(CuCl₂) | Cu 寄与 |
|---|---|---|---|---|---|---|---|
| IdealTemkin | 0.984 | 7.00e-3 | 4.28e-1 | 1.69e-14 | 2.40e-11 | **1.55e-6** | CuCl₂ 100% |
| RegularSolution | 0.915 | 1.21e-2 | 7.38e-1 | 2.91e-14 | 1.23e-10 | **2.68e-6** | CuCl₂ 100% |

CuCl(g) と Cu₃Cl₃(g) だけなら寿命は 1.8×10⁷ 〜 9.2×10⁷ h。
**寿命制約は事実上 CuCl₂(g) だけへの制約になっている。**

RegularSolution は `{('CuCl','KCl'): -24400, ('CuCl2','KCl'): +16700, ('CuCl','CuCl2'): -6800}` [J/mol]。
W(CuCl₂,KCl) が正なので KCl は CuCl₂ を**不安定化**し、理想 Temkin より a(CuCl₂) が
大きくなる。したがって RegularSolution 側のほうが制約は厳しい（下限が高い）。

### 結果

| 活量モデル | ΔHf₂₉₈(CuCl₂,g) 下限 [kJ/mol] |
|---|---|
| IdealTemkin | **≥ -38.8** |
| RegularSolution | **≥ -36.1** |
| （参考）旧 S/Cp のまま IdealTemkin | ≥ -44.1 |
| （参考）旧 S/Cp のまま RegularSolution | ≥ -41.3 |

**幅として ΔHf₂₉₈ ≥ -38.8 〜 -36.1 kJ/mol、丸めて ≥ -39 kJ/mol。**

現行暫定値 **-43.0 はこの下限を 4.2 〜 6.9 kJ/mol 割っている**。
新しい S°/Cp と組み合わせると、参照触媒の寿命は 2.5×10³ 〜 4.3×10³ h と
予測され、実測の 9600 h に 2.3 〜 3.9 倍足りない。

なお旧 S/Cp のままだと -43.0 は下限（-44.1）のすぐ内側にあり、
「暫定値が偶然もっともらしく見えていた」のは S° を 7.5 J/mol/K 低く取っていた
ためだと分かる。

### 下限のロバストネス

寿命目標を振ったとき（新 S/Cp）:

| 目標 t_half | IdealTemkin | RegularSolution |
|---|---|---|
| 4800 h | -42.4 | -39.6 |
| 9600 h | -38.8 | -36.1 |
| 19200 h | -35.3 | -32.5 |

寿命 2 倍 ⇔ 下限 3.5 kJ/mol。GHSV や Cu 装填量の想定も同じ係数で効く
（GHSV が半分なら下限は 3.5 kJ/mol 下がる）。

S°₂₉₈ を ±5 J/mol/K 振ったとき:

| S°₂₉₈ | IdealTemkin | RegularSolution |
|---|---|---|
| 280.5 | -41.9 | -39.1 |
| 285.5 | -38.8 | -36.1 |
| 290.5 | -35.8 | -33.0 |

### この下限の性格（重要）

これは **ΔHf 単独の熱力学的下限ではなく、「ΔHf + 活量モデル + 反応器条件」の
連立に対する制約**である。等価な逃げ道が 3 つある:

1. ΔHf(CuCl₂,g) が実際に -39 kJ/mol より正（＝ CuCl₂(g) がもっと不安定）
2. 融液中の a(CuCl₂) が両モデルより 3-4 倍小さい（クロロ銅酸錯体 CuCl₃²⁻/CuCl₄²⁻
   による安定化。CLAUDE.md の未決着事項そのもの）
3. Cu が担体（γ-Al₂O₃）や K-Cu 複塩に固定されていて融液の Cu 活量が下がっている

**したがって「ΔHf ≥ -39 kJ/mol」は断定ではなく、「-43.0 を使い続けるなら
a(CuCl₂) を 3-4 倍下げる機構を同時に持ち込まないと実測寿命と矛盾する」
という形の拘束**として読むのが正しい。逆に言えば、ΔHf の一次資料
（Wächter & Schäfer 1980 など）が取れれば、この不等式は
**活量モデルの較正データ**に転用できる。

---

## 5. 感度と未取得

### S°₂₉₈ の感度（298.15 K、単独で振った場合）

| 振った量 | 変化 | ΔS°₂₉₈ [J/mol/K] |
|---|---|---|
| ν₂ 95.81 → 127 cm⁻¹ | 変角が硬い側 | **-4.58** |
| ν₂ 95.81 → 80 cm⁻¹ | 変角が柔らかい側 | **+2.95** |
| ²Σ_g⁺ 1911 → 99 cm⁻¹（ab initio 説） | | **+4.53** |
| ²Σ_g⁺ 1911 → 500 cm⁻¹ | | +2.03 |
| ²Σ_g⁺ 1911 → 1600 / 2400 cm⁻¹ | | +0.02 / -0.01 |
| スピン軌道分裂を無視（g = 4, E = 0） | | +3.27 |
| r 2.0352 → 2.10 Å | | +0.52 |

**支配的な不確かさは (a) ν₂ と (b) ²Σ_g⁺ の位置**。r は 3% 動かしても 0.5 しか効かない。
両者を合わせて **S°₂₉₈ = 285.5 ± 4 J/mol/K** とする（653 K では ±4 が p に ×0.79-1.27）。

### 未取得・要注意

1. **ν₂ が調和振動数 ω₂ か観測振動間隔かが確認できていない。** Beattie 1997 は
   「ω₁, ω₂, ω₃ と回転・非調和定数を報告」とあるが、本文が有料で確認できず、
   NIST WebBook の表の 95.81 がどちらかは不明。加えて X ²Π_g では
   **Renner-Teller 相互作用**により変角準位が単純な調和梯子にならない
   （NIST WebBook の (0 2l 0) 準位に関する記載、Renner-Teller と K-doubling）。
   ν₂ が 96 cm⁻¹ と小さいので、この効果は S に数 J/mol/K の系統誤差を残す。
2. **²Σ_g⁺ の気相値が未測定。** Ne 1910.9 / Ar 1616.5 しかない。マトリックスシフトは
   ²Π_1/2 で見ると Ne -8.8 / Ar -179 なので Ne は信頼できるが、
   ab initio が ~100 cm⁻¹ を出している以上、割り切れていない。
   → もし ²Σ_g⁺ が実際に低ければ S はさらに +4.5、下限はさらに +2.8 kJ/mol 上がる
   （つまり本ノートの下限は**保守側**）。
3. **ΔHf₂₉₈ の一次資料は今回取れていない。** JANAF / NIST WebBook に気相
   CuCl₂ の熱化学値は無い。OA 範囲では Wächter & Schäfer (1980) も
   Gurvich/IVTANTHERMO も参照できなかった。ΔHf は EST のまま残す。
4. **非調和性・振動回転相互作用**は入れていない。伸縮モードでは 653 K でも
   0.2 J/mol/K 程度の効果。無視してよい。
5. **Cu₂Cl₄(g) 等の Cu(II) 二量体**は DB に無い。CuCl 系で三量体が支配的だった
   ことを考えると、Cu(II) 側にも会合種がある可能性は残る（未検討）。

---

## 6. `data.py` への反映提案

S°/Cp は分光定数から一意に決まるので **計算値で置換**、ΔHf は
一次資料が無いので **EST のまま、下限を注記**する形を提案する。

```python
DB.add(
    Species(
        "CuCl2(g)",
        "g",
        -43.0,          # ★ EST のまま。docs/cucl2g_statmech.md の下限 -39 を割っている
        285.5,          # 統計力学 (剛体回転子+調和振動子+電子項)
        (64.845, 0.3103, -3.8724),   # 同 MK フィット (298-1000 K, 最大残差 0.06)
        {"Cu": 1, "Cl": 2},
        EST,
        "S298/Cp は分光定数からの統計力学計算 (2026-08-10, docs/cucl2g_statmech.md)。"
        "r0=203.52 pm (Beattie 1997, Inorg. Chem. 36, 3207), "
        "nu=(369.36, 95.81x2, 522.18) cm-1, 電子項 2Pi3/2(g=2,0) / 2Pi1/2(g=2,482.9) / "
        "2Sg+(g=2,1911) / 2D5/2(g=2,6877) / 2D3/2(g=2,8700) [NIST WebBook B1001128, "
        "Hodges & Ross JCP 127 (2007) 024309]。同コードで CuCl(g) S298 を "
        "236.96 (JANAF 237.20) と再現。S の不確かさ ±4 (支配は nu2 と 2Sg+ の位置)。"
        "★ dHf298 は依然 EST: JANAF/WebBook に気相データ無し。参照触媒の 9600 h 寿命は "
        "dHf298 >= -38.8 (IdealTemkin) 〜 -36.1 (RegularSolution) kJ/mol を要求する。"
        "-43.0 を使い続けるなら a(CuCl2) を 3-4 倍下げる機構 (クロロ銅酸錯体) が要る。"
        "Waechter & Schaefer (1980) 要取得",
    )
)
```

**手順の注意（CLAUDE.md のルール）**: この差し替えは揮発の結論を約 3 倍動かすので、
**先にコミットしてから図を再生成すること**。既存の図はすべて旧 S/Cp に対応している。

あわせて確認すべきテスト:
- `tests/test_model.py::test_ideal_model_reproduces_observed_lifetime`（xfail）は
  この変更で**さらに外れる方向**に動く（予測寿命 1.19e4 h → 4.25e3 h）。
  xfail の理由書きを「Cu₃Cl₃ データ誤差」から「a(CuCl₂) の錯体安定化 + ΔHf(CuCl₂,g)」に
  更新するのが妥当。

---

## 付録: 再現コード（要点のみ）

```python
import numpy as np
R, NA, KB, H, C, AMU = 8.314462618, 6.02214076e23, 1.380649e-23, 6.62607015e-34, 2.99792458e10, 1.66053906660e-27
CM_TO_K = H * C / KB

def s_trans(M, T, P=101325.0):
    m = M * 1e-3 / NA
    q = (2*np.pi*m*KB*T/H**2)**1.5 * KB*T/P
    return R * (np.log(q) + 2.5)

def s_rot_lin(I, T, sigma):
    return R * (np.log(8*np.pi**2*I*KB*T/(sigma*H**2)) + 1.0)

def vib(nus, T):            # nus = [(cm-1, degeneracy), ...]
    S = Cv = 0.0
    for nu, g in nus:
        x = nu * CM_TO_K / T
        S  += g*R*(x/np.expm1(x) - np.log(-np.expm1(-x)))
        Cv += g*R*x**2*np.exp(x)/np.expm1(x)**2
    return S, Cv

def elec(levels, T):        # levels = [(cm-1, g), ...]
    b = np.array([e for e, _ in levels])*CM_TO_K/T
    g = np.array([gg for _, gg in levels], float)
    w = g*np.exp(-b); q = w.sum()
    m1 = (w*b).sum()/q; m2 = (w*b**2).sum()/q
    return R*(np.log(q)+m1), R*(m2-m1**2)

M_CU, M_CL, r = 63.546, 35.453, 2.0352          # Å
I = 2 * M_CL*AMU * (r*1e-10)**2                  # 直線 XY2
NUS    = [(369.36, 1), (95.81, 2), (522.18, 1)]
LEVELS = [(0.0, 2), (482.9, 2), (1910.9, 2), (6877.0, 2), (8700.0, 2)]

def cucl2(T):
    Sv, Cvv = vib(NUS, T); Se, Cve = elec(LEVELS, T)
    S  = s_trans(M_CU+2*M_CL, T) + s_rot_lin(I, T, 2) + Sv + Se
    Cp = 3.5*R + Cvv + Cve                       # 直線分子: 3/2R + R + R
    return S, Cp

# 検証: CuCl(g)  ->  S298 = 236.96 (JANAF 237.20)
I1 = (M_CU*M_CL/(M_CU+M_CL))*AMU*(2.0512e-10)**2
S_cucl = s_trans(M_CU+M_CL, 298.15) + s_rot_lin(I1, 298.15, 1) + vib([(415.3,1)], 298.15)[0]
```

`ΔHf` 下限の逆算（リポジトリのコードを使う）:

```python
from scipy.optimize import brentq
from deacon_thermo.data import DB
from deacon_thermo.melt import IdealTemkin, RegularSolution, redox_split
from deacon_thermo.species import Confidence, Species
from deacon_thermo.volatility import lifetime

T, P_CL2, DIL, TARGET = 613.15, 0.32, {"KCl": 1.0, "LaCl3": 0.3}, 9600.0

def bound(model, S=285.46, cp=(64.845, 0.3103, -3.8724)):
    def f(dHf):
        DB["CuCl2(g)"] = Species("CuCl2(g)", "g", dHf, S, cp, {"Cu": 1, "Cl": 2},
                                 Confidence.ESTIMATE, "scan")
        _, melt = redox_split(1.0, P_CL2, T, DIL, model)
        return np.log(lifetime(melt, T)[0]) - np.log(TARGET)
    return brentq(f, -300.0, 600.0, xtol=1e-3)

bound(IdealTemkin())                                            # -> -38.8
bound(RegularSolution({("CuCl","KCl"): -24400.0,
                       ("CuCl2","KCl"): 16700.0,
                       ("CuCl","CuCl2"): -6800.0}))             # -> -36.1
```

## 参照 URL 一覧（すべて 2026-08-10 取得）

- NIST Chemistry WebBook, CuCl₂ 振動・電子エネルギー準位: https://webbook.nist.gov/cgi/cbook.cgi?ID=B1001128&Mask=800
- NIST CCCBDB, CuCl 実験データ: https://cccbdb.nist.gov/exp2x.asp?casno=7758896
- Beattie et al., *Inorg. Chem.* 36 (1997) 3207（書誌・要旨）: https://pubs.acs.org/doi/10.1021/ic961509r , https://pubmed.ncbi.nlm.nih.gov/11669980/
- Hodges & Ross, *J. Chem. Phys.* 127 (2007) 024309（要旨）: https://pubs.aip.org/aip/jcp/article/127/2/024309/919449 , https://www.ncbi.nlm.nih.gov/pubmed/17640130
- Ab initio study on the spectroscopy of CuCl₂ I / II（要旨）: https://pubmed.ncbi.nlm.nih.gov/15268475/ , https://pubmed.ncbi.nlm.nih.gov/15945683/
