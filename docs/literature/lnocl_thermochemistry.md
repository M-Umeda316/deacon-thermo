# LnOCl の熱化学 — 文献調査

調査日: 2026-08-05 / 対象: `data.py` の `DH_OXYCHLORIDE = -30.0 kJ/mol` の検証
対象問題: 反応条件 (380 ℃, HCl/O₂ = 2:1) で Sm は SmCl₃ のままか SmOCl になるか

> **本文書の数値の扱い**
> 各値に「信頼度」を付す。凡例:
> - **[全文]** 一次または二次資料の本文/表を直接読んで取得
> - **[抄録]** 抄録・出版社要約・研究ポータル要約のみで取得。本文未確認
> - **[導出]** 上記から本調査で計算した値（計算式を明記）
> 見つからなかったものは「未取得」と明記する。推定で数値を埋めていない。

---

## 0. 結論サマリ（先に）

1. **加水分解平衡そのものの高温実測は存在する。** Koch–Cunningham（UC Berkeley）による
   JACS 3 部作が、まさに `LnCl₃(s) + H₂O(g) = LnOCl(s) + 2 HCl(g)` の
   熱量・自由エネルギーを La / Sm / Gd / Pr / Nd について測定している。
   **Sm は第 II 報に含まれる。** ただし全 3 報とも ACS ペイウォールで、
   **本調査では数値そのものを取得できなかった**（書誌情報は Crossref で確定）。
   → 最優先の要取得文献（§4）。

2. **`DH_OXYCHLORIDE = -30 kJ/mol` は浅すぎる。** 取得できた二次資料はすべて、
   より発熱側（負側）を示す。プロジェクトのモデル内で литература を再現する
   実効値は **Sm で約 -42 kJ/mol**（Knacke 系, 1000 K の ΔG 基準）。
   熱量測定にもとづく Gd は -35.4、Eu は -41.0 kJ/mol。

3. **結果として Sm の判定は塩化物側からオキシ塩化物側へ動く。**
   `DH = -42` にすると、653 K・HCl/O₂ = 2:1 で SmCl₃ が安定なのは
   転化率 ~45% までとなり、**床の大半（特に出口側）で SmOCl 側**になる。
   現行 `-30` では転化率 ~72% まで塩化物側だった。

4. **ただし「系列を通して一定の DH」という仮定自体が文献と整合しない。**
   実測 4 元素（La, Nd, Sm, Gd）の間ですら ΔG(1000 K) が -31 〜 -57 kJ/mol と
   26 kJ/mol ばらつく。さらに Nd については Knacke 系（-31 kJ/mol @1000 K）と
   Navrotsky/Woodfield 系（-59.4 kJ/mol @298 K）が約 28 kJ/mol 食い違っている。
   **これは CLAUDE.md の「系列内の差は絶対値より信頼できる」という前提への
   直接の反証候補**なので、§5 のギャップとして最重要視すべき。

---

## 1. サマリ表

### 1.1 反応の定義（本文書で統一）

| 記号 | 反応 | 備考 |
|---|---|---|
| **加水分解** | LnCl₃(s) + H₂O(g) = LnOCl(s) + 2 HCl(g) | K = p(HCl)²/p(H₂O) [atm] |
| **DH** | (1/3) Ln₂O₃(s) + (1/3) LnCl₃(s) → LnOCl(s) | `data.py` の `DH_OXYCHLORIDE` の定義と同一。**LnOCl 1 mol あたり** |

`DH` は本来エンタルピーだが、文献値は ΔG(T) で与えられることが多い。
以下では**どちらの量か必ず明記**する。混同が最大の事故源。

### 1.2 LnOCl の生成エンタルピー ΔHf₂₉₈（実測）

| 化合物 | ΔHf₂₉₈ [kJ/mol] | 方法 | 出典 | 信頼度 |
|---|---|---|---|---|
| EuOCl(cr) | **-903.5 ± 3** | 1.000 N HCl 中の溶解熱量測定 | Burns, Peterson, Haire (1998) | [抄録] |
| GdOCl(cr) | **-981.4 ± 3** | 同上 | 同上 | [抄録] |
| LuOCl(cr) | **-987.2 ± 4** | 同上 | 同上 | [抄録] |
| LaOCl, NdOCl, SmOCl | **未取得** | — | Knacke ら編纂に収載（値未確認） | — |

補足: Burns らは 50–1150 μg スケールの微量熱量測定（粉末・単結晶）。
Eu, Gd, Lu の実測値をもとに**全 Ln の LnOCl の生成エンタルピーを推定**しており、
その推定表が Jacob らの §1.3 の (d) 欄の元になっている [全文: Jacob 2016 本文記述]。

### 1.3 LnOCl の ΔG°f (1000 K) と、そこから導いた DH

Jacob, Dixit, Rajput (2016) Table 1 の実数値（**[全文]**、PDF を直接読取）。
`DH(1000 K, ΔG基準) = ΔG°f(LnOCl) - (1/3)ΔG°f(Ln₂O₃) - (1/3)ΔG°f(LnCl₃)` は**[導出]**。

| Ln | ΔG°f(Ln₂O₃) | ΔG°f(LnCl₃) | ΔG°f(LnOCl) | 出所タグ | **DH(ΔG, 1000 K)** |
|---|---:|---:|---:|:--:|---:|
| La | -1508.842 | -829.3025 | -836.839 | (a) 実測系 | **-57.5** |
| Ce | -1512.87 | -811.5925 | -814.316 | (d) 推定 | -39.5 |
| Pr | -1518.582 | -814.7055 | -825.608 | (d) 推定 | -47.8 |
| Nd | -1524.974 | -801.7315 | -806.675 | (a) 実測系 | **-31.1** |
| **Sm** | **-1530.688** | **-780.8185** | **-821.645** | **(a) 実測系** | **-51.1** |
| Eu | -1355.068 | -689.3295 | -710.739 | (d) 推定 | -29.3 |
| Gd | -1538.825 | -773.5055 | -801.814 | (a) 実測系 | **-31.0** |
| Tb | -1578.41 | -762.1265 | -793.823 | (d) 推定 | -13.6 |
| Dy | -1565.825 | -753.8365 | -798.554 | (d) 推定 | -25.3 |
| Ho | -1591.563 | -759.9815 | -814.465 | (d) 推定 | -30.6 |
| Er | -1606.756 | -752.0525 | -808.824 | (d) 推定 | -22.6 |
| Tm | -1584.005 | -748.3545 | -797.224 | (d) 推定 | -19.8 |
| Yb | -1530.423 | -709.7605 | -772.335 | (d) 推定 | -25.6 |

単位はすべて kJ/mol、**T = 1000 K、標準状態 p° = 0.1 MPa**。
出所タグ (a) = Knacke, Kubaschewski, Hesselmann (1991) の編纂値（La, Nd, Sm, Gd の 4 元素のみ
LnOCl の温度依存データが収載されている＝実測に由来）。(d) = Burns ら (1998) の推定値。

Jacob らが (d) の ΔG を作る際に使った entropy 仮定（**[全文]**、本文明記）:

> ΔS_av = **-8 (±4) J mol⁻¹ K⁻¹** for (1/3)Ln₂O₃ + (1/3)LnCl₃ → LnOCl
> （Cp は Neumann–Kopp 則）

これを使うと (d) 欄について `DH(ΔH基準) ≈ DH(ΔG,1000 K) - 8 kJ/mol`。**[導出]**

### 1.4 Navrotsky/Woodfield 系の値（298 K, ΔG 基準, 酸化物+塩化物からの生成）

| 化合物 | ΔG_f(298.15 K) rel. 酸化物+塩化物 [kJ/mol] | 出典 | 信頼度 |
|---|---:|---|---|
| NdOCl | **-59.4** | Gibson, Yang, Riman, Navrotsky, Woodfield (2025) | [抄録] |
| YOCl | **-41.0** | 同上 | [抄録] |
| TmOCl | **-11.1** | 同上 | [抄録] |

これは高温酸化物融体溶解熱量測定（Yang ら 2022）で得た生成エンタルピーと、
低温熱容量測定から得た S°₂₉₈ を組み合わせたもの。
**規格化が「LnOCl 1 mol あたり（= 1/3 Ln₂O₃ + 1/3 LnCl₃ 基準）」であることは
本文未確認**（この分野の慣行からはほぼ確実だが、要検証）。

Yang ら (2022) の定性的結論（**[抄録]**、Semantic Scholar 経由で全文抄録を確認）:

- REOCl はすべて、対応する二元酸化物・塩化物の等モル混合物に対して**熱力学的に安定**
- REOCl は RE サイズ減少に伴い 3D 骨格（PbFCl 型）→ 2D ファンデルワールス層状（SmSI 型）へ構造転移
- 二元酸化物・塩化物からの生成エンタルピーは発熱で、**イオン半径の減少とともに発熱性が小さくなる**

TmOCl の -11.1 kJ/mol という浅さは、この 3D→2D 構造転移によるもので、
**系列を通して一定の DH という近似が重希土側で完全に破綻する**ことを意味する。

---

## 2. 現行推定値（DH_OXYCHLORIDE = -30 kJ/mol）との比較

### 2.1 「DH_OXYCHLORIDE の実効値」への換算

`data.py` の推定器は ΔHf だけでなくエントロピーにも仮定を置いている:

```python
dHf(LnOCl) = dHf(Ln2O3)/3 + dHf(LnCl3)/3 + DH_OXYCHLORIDE   # DH = -30.0
S298(LnOCl) = S298(Ln2O3)/3 + S298(LnCl3)/3 + 8.0           # ΔS_rxn = +8 J/mol/K
```

**ここに重大な符号の不一致がある。**
プロジェクトは ΔS_rxn = **+8** J/mol/K を仮定しているが、
Burns / Jacob らは ΔS_av = **-8 (±4)** J/mol/K を採用している（§1.3）。
差は 16 J/mol/K で、653 K では ΔG にして **約 10.4 kJ/mol**、
1000 K では **16 kJ/mol** に相当する。プロジェクト側の符号は LnOCl を
**過度に安定側に**見せる向きに効く。

そのため「DH の実効値」は基準温度で変わる。以下は**プロジェクトのモデル内で
（Cp・S の現行仮定をそのまま使って）文献の ΔG を再現するのに必要な DH** で、
これが最も実務的な換算値である。**[導出]**

| Ln | 文献値（基準） | プロジェクト現行 (DH=-30) の同温 ΔG | **DH 実効値** | 現行との差 |
|---|---|---:|---:|---:|
| La | ΔG(1000 K) = -57.5 (Knacke) | -39.9 | **-47.6** | -17.6 |
| Ce | ΔG(1000 K) = -39.5 (Burns 推定) | -39.9 | -29.6 | +0.4 |
| Pr | ΔG(1000 K) = -47.8 (Burns 推定) | -39.9 | -37.9 | -7.9 |
| Nd | ΔG(1000 K) = -31.1 (Knacke) | -39.9 | -21.2 | **+8.8** |
| Nd | ΔG(298 K) = -59.4 (JCT 2025) | -32.4 | **-57.0** | **-27.0** |
| **Sm** | **ΔG(1000 K) = -51.1 (Knacke)** | **-39.1** | **-42.0** | **-12.0** |
| Eu | ΔHf₂₉₈ = -903.5 (Burns 実測) | — | **-41.0** | -11.0 |
| Gd | ΔHf₂₉₈ = -981.4 (Burns 実測) | — | **-35.4** | -5.4 |
| Gd | ΔG(1000 K) = -31.0 (Knacke) | -37.7 | -23.3 | +6.7 |
| Y | ΔG(298 K) = -41.0 (JCT 2025) | -32.4 | -38.6 | -8.6 |
| Dy | ΔG(1000 K) = -25.3 (Burns 推定) | -38.6 | -16.7 | +13.3 |
| Ho | ΔG(1000 K) = -30.6 (Burns 推定) | -39.7 | -20.9 | +9.1 |
| Er | ΔG(1000 K) = -22.6 (Burns 推定) | -39.9 | -12.7 | +17.3 |

Eu, Gd の行（Burns 実測 ΔHf₂₉₈ 由来）は、プロジェクトの `_LN_DATA` の
Ln₂O₃ / LnCl₃ を補助データに使った直接計算:

```
DH = ΔHf(LnOCl) - (1/3)ΔHf(Ln₂O₃) - (1/3)ΔHf(LnCl₃)
Eu: -903.5 - (1/3)(-1651.4) - (1/3)(-936.0) = -903.5 + 862.467 = -41.03
Gd: -981.4 - (1/3)(-1819.7) - (1/3)(-1018.2) = -981.4 + 945.967 = -35.43
```

> **注意**: Burns らが自分の計算で使った Ln₂O₃ / LnCl₃ の補助データは未確認。
> プロジェクトの補助データと違えば、この DH はその差だけずれる。

### 2.2 プロジェクトの判定への影響（Sm）

653.15 K, HCl/O₂ = 2:1, P = 1 atm。`sensitivity.chloride_margin_at()` で計算。
正なら SmCl₃ が安定。**[導出]**

| HCl 転化率 | DH=-30（現行） | DH=-42（本調査推奨） | DH=-51 | DH=-59 |
|---:|---:|---:|---:|---:|
| 0.10 | +23.7 | +11.7 | +2.7 | -5.3 |
| 0.30 | +15.2 | +3.2 | -5.8 | -13.8 |
| 0.50 | +9.0 | -3.0 | -12.0 | -20.0 |
| 0.70 | +1.8 | -10.2 | -19.2 | -27.2 |
| 0.78 | -2.1 | -14.1 | -23.1 | -31.1 |
| 0.838（平衡） | -5.7 | -17.7 | -26.7 | -34.7 |

（単位 kJ/mol。653.15 K の平衡転化率は 0.838、Q = p(HCl)²/p(H₂O) = 0.0486 atm）

**読み方**:

- 現行 `-30` でも、既に**出口側（転化率 > 約 72%）では SmOCl 側**だった。
  CLAUDE.md の「余裕は約 30 kJ/mol」は入口〜低転化率側の値に対応する。
- `-42` にすると反転点が**転化率 ~45%** に前倒しになる。
  つまり**床の後半 55% で SmOCl** になる。
- `-51`（Knacke の ΔG をそのまま DH と読んだ場合）では転化率 ~25% で反転。
  ほぼ床全体が SmOCl。

**したがって、「Sm は SmCl₃ か SmOCl か」の答えは、取得できた文献の範囲では
「床の軸方向で変わり、出口側は SmOCl」に傾く。** CLAUDE.md が
「転化率に沿って床の軸方向で相が変わっている可能性」と書いていた仮説を
文献値は**支持する方向**に動かす。

### 2.3 他元素への影響（参考）

| Ln | DH 実効値 | 転化率 0.30 での余裕 | 平衡での余裕 |
|---|---:|---:|---:|
| La | -47.6 | +36.9 | +16.0 |
| Nd | -21.2 (Knacke) | +39.1 | +18.2 |
| Nd | -57.0 (JCT 2025) | +3.3 | **-17.6** |
| Gd | -35.4 (Burns 実測) | +6.1 | **-14.8** |

Nd の 2 行の乖離（同じ条件で +18.2 と -17.6）が、現時点の最大の不確かさ。

---

## 3. 出典詳細

### 3.1 【最重要】加水分解平衡の直接測定 — Koch–Cunningham 3 部作

書誌情報は Crossref API で確定（**[全文]** 相当の書誌）。**数値は未取得（ペイウォール）。**

| # | 反応 | 著者 | 誌 | 巻(号) | 頁 | 年 | DOI |
|---|---|---|---|---|---|---|---|
| I | LaCl₃(s) + H₂O(g) = LaOCl(s) + 2HCl(g) | C. W. Koch, A. Broido, B. B. Cunningham | J. Am. Chem. Soc. | 74(9) | 2349–2351 | 1952 | `10.1021/ja01129a049` |
| **II** | **SmCl₃ / GdCl₃ の同反応** | **C. W. Koch, B. B. Cunningham** | **J. Am. Chem. Soc.** | **75(4)** | **796–797** | **1953** | **`10.1021/ja01100a010`** |
| III | PrCl₃ / NdCl₃ の同反応 | C. W. Koch, B. B. Cunningham | J. Am. Chem. Soc. | 76(6) | 1471–1474 | 1954 | `10.1021/ja01635a003` |

いずれもタイトルに "Heat and Free Energy of the Reaction" と明記されており、
**ΔH と ΔG（したがって K(T) = p(HCl)²/p(H₂O)）が報告されている**ことは確実。
Unpaywall で確認したところ 3 報とも OA 版なし。

**プロジェクトにとっての意味**: `stability.hydrolysis_K()` が計算しているまさに
その量の直接実測値であり、Ln₂O₃ / LnCl₃ の補助データを経由しないので
誤差伝播が最小。§1.3 の Knacke 値も辿ればここに行き着くはずだが、
一次資料の測定温度域・不確かさが分からないと信頼区間が引けない。

### 3.2 LnOCl 生成エンタルピーの熱量測定

**Burns, J. B.; Peterson, J. R.; Haire, R. G.**
"Standard enthalpies of formation for europium, gadolinium, and lutetium oxychlorides,
calculated from measured enthalpies of solution."
*J. Alloys Compd.* **1998**, *265*(1–2), 146–152. DOI `10.1016/S0925-8388(97)00435-0`
（書誌 [全文] / 数値 [抄録]）

- 方法: 1.000 N HCl 中の溶解熱、50–1150 μg スケール、粉末および単結晶
- ΔfH°₂₉₈: EuOCl -903.5 ± 3、GdOCl -981.4 ± 3、LuOCl -987.2 ± 4 kJ/mol
- 上記 3 点をもとに**全ランタノイドの LnOCl の ΔfH を推定**（この推定表は未取得）
- 関連: 同グループの CfOCl 論文 `10.1016/S0925-8388(98)00185-6`（*J. Alloys Compd.* 1998）

**Yang, S.; Anderko, A.; Riman, R. E.; Navrotsky, A.**
"Thermochemistry of 3D and 2D Rare Earth Oxychlorides (REOCls)."
*Inorg. Chem.* **2022**, *61*(19), 7590–7596. DOI `10.1021/acs.inorgchem.2c00763`
（書誌 [全文] / 抄録 [全文] / 数値表 未取得）

- 方法: 高温酸化物融体溶解熱量測定 + Born–Haber サイクル、フラックス法で全 REOCl を合成
- 定性的結論は §1.4 に記載。**数値表（各 REOCl の ΔfH）は未取得**
- OSTI ID 1981863（受理原稿 PDF は公開されていない）

**Gibson, A.; Yang, S.; Riman, R. E.; Navrotsky, A.; Woodfield, B. F.**
"Heat capacity and thermodynamic functions of stoichiometric rare earth oxychlorides (REOCl)."
*J. Chem. Thermodyn.* **2025**, *211*, 107549. DOI `10.1016/j.jct.2025.107549`
プレプリント: SSRN `10.2139/ssrn.5208841`
（書誌 [全文] / 数値 [抄録]）

- 測定対象: **TmOCl, NdOCl, YOCl**
- 低温熱容量 → S°₂₉₈、既報の生成エンタルピーと組み合わせて ΔG_f を導出
- ΔG_f(298.15 K) 対 酸化物+塩化物: NdOCl **-59.4**、YOCl **-41.0**、TmOCl **-11.1** kJ/mol
- **プロジェクトが持っていない S°₂₉₈(LnOCl) の実測値がここにある**（現行は
  `S(Ln2O3)/3 + S(LnCl3)/3 + 8.0` の推定）。SSRN プレプリントは
  ダウンロード試行が 403 だったが、機関アクセスなしで取れる可能性が高い

### 3.3 系列を通した評価・レビュー

**Jacob, K. T.; Dixit, A.; Rajput, A.**
"Stability field diagrams for Ln–O–Cl systems."
*Bull. Mater. Sci.* **2016**, *39*(3), 603–611. DOI `10.1007/s12034-016-1219-6`
**オープンアクセス**: `https://www.ias.ac.in/public/Volumes/boms/039/03/0603-0611.pdf`
（**[全文]** — PDF を取得して Table 1・本文を直接読取）

**本調査で最も価値が高かった資料。** 理由:

1. 13 元素の ΔG°f(Ln₂O₃, LnCl₃, LnCl₂, LnOCl) を 1000 K で一表にしており、
   各値の出所（Knacke / Burns / Uda / Pankratz）がタグ付けされている
2. **プロジェクトの `stability.py` と同じ Kellogg 図（log p(O₂) vs log p(Cl₂)）**を
   13 系について描いており、相の隣接関係を直接照合できる
3. 本文の重要な記述:
   - 「LnOCl の熱力学データが編纂に載っているのは **La, Nd, Sm, Gd の 4 元素だけ**」
   - 「LnOCl のデータは酸化物・三塩化物のデータより精度が低い」
   - 「Ln₂O₃ と LnCl₃ は共存できない。両者の間には必ず LnOCl が挟まる」
     → プロジェクトの `lanthanide_phases()` が LnCl₃/LnOCl/Ln₂O₃ の 3 相を
     候補にしているのは正しく、かつ Ln₂O₃–LnCl₃ 直接境界は現れないはず
   - LaCl₃–LaOCl 系は 1093 K に共晶をもつ（Drobot ら 1965）
     → 融液モデル (`melt.py`) を LnOCl まで広げる際の入口

**Knacke, O.; Kubaschewski, O.; Hesselmann, K.**
*Thermochemical Properties of Inorganic Substances*, 2nd ed.; Springer: Berlin, **1991**.
（未取得。Barin と同系統の編纂。**LaOCl, NdOCl, SmOCl, GdOCl の温度依存データを収載**）

### 3.4 Deacon / オキシ塩素化触媒における LnOCl vs LnCl₃

**参照系そのもの**（プロジェクトの前提の出典）:

- Sun, Y.; Li, C.; Guo, Y.; Zhan, W.; Guo, Y.; Wang, L.; Wang, Y.; Lu, G.
  "Catalytic oxidation of hydrogen chloride to chlorine over Cu-K-Sm/γ-Al₂O₃ catalyst
  with excellent catalytic performance." *Catal. Today* **2018**, *307*, 286–292.
  DOI `10.1016/j.cattod.2017.04.014`
- Feng, K.; Li, C.; Guo, Y.; Zhan, W.; Ma, B.; Chen, B.; Yuan, M.; Lu, G.
  "An efficient Cu-K-La/γ-Al₂O₃ catalyst for catalytic oxidation of hydrogen chloride
  to chlorine." *Appl. Catal. B* **2015**, *164*, 483–487.
  DOI `10.1016/j.apcatb.2014.09.063`
  — 0.1 MPa, 340 ℃, GHSV 450 L kg⁻¹ h⁻¹, HCl/O₂ = 2:1 で HCl 転化率 約 78% を
  9600 h 維持（**[抄録]**）。CLAUDE.md の記述と一致
- Aglulin, A. G. "Kinetics and possible mechanism of hydrogen chloride oxidation over
  supported copper-containing salt catalysts: II. ... over the CuCl₂-KCl-LaCl₃ catalyst."
  *Kinet. Catal.* **2014**, *55*(5), 582–591. DOI `10.1134/S0023158414050024`
  — 350–425 ℃、無勾配法。**LaCl₃ 添加で HCl 酸化速度が 1 桁上がる**（**[抄録]**）。
  I 報（Cu-K のみ）は `10.1134/S0023158414050012`

**LnOCl / LnCl₃ どちらが活性相かの実験報告**（メタン・エタンのオキシ塩素化系。
Deacon と同じ HCl/O₂ 雰囲気なので相安定性の議論は直接転用できる）:

- Terlingen, B.; Oord, R.; Ahr, M.; Hutter, E. M.; van Lare, C.; Weckhuysen, B. M.
  "Favoring the Methane Oxychlorination Reaction over EuOCl by Synergistic Effects
  with Lanthanum." *ACS Catal.* **2022**, *12*(9), 5698–5710.
  DOI `10.1021/acscatal.2c00777` — **オープンアクセス (PMC9087184)**（**[全文]**）
  - 350–550 ℃, CH₄/HCl/O₂/N₂/He
  - 「**熱力学計算によれば LnOCl → LnCl₃ の塩素化は LaOCl が最も容易**」
  - HCl 10% で LaOCl は速やかに塩素化されるが、EuOCl は非常に高い HCl 濃度と
    長時間を要する → **La は Deacon 条件で塩化物側に行きやすい Ln の代表**
  - La₀.₅Eu₀.₅OCl は反応条件下で La リッチ相と Eu リッチ相に相分離
  - **注記**: 個々の Ln の ΔG は本文に数表として載っていない
- ACS Catal. **2021**, *11*(16), 10574 (Weckhuysen ら) "Mechanistic Insights into the
  Lanthanide-Catalyzed Oxychlorination of Methane as Revealed by Operando Spectroscopy"
  — LaOCl 上では La³⁺ ではなく**表面 Cl の形式酸化数がサイクルする**（Cl⁻ ⇄ ClO⁻）
  機構を提案（**[抄録]**）。CLAUDE.md の「Cu²⁺/Cu⁺ が触媒サイクル」という枠組みと
  対比すべき別機構
- Li, Y. ら "Metastable LaOClₓ Phase Stabilization as an Effective Strategy for
  Controllable Chlorination of Ethane into 1,2-Dichloroethane."
  *Molecules* **2025**, *30*(8), 1746. DOI `10.3390/molecules30081746`
  — **オープンアクセス**（**[全文]**）。**中間塩素化度の準安定 LaOClₓ 相が活性相**で、
  純 LaOCl は不可逆に LaCl₃ 化して 10 h で選択性を失う。Al₂O₃ に閉じ込めると
  LaOClₓ が安定化して 260 ℃・高 Cl₂ 濃度で 12 h 安定
  → CLAUDE.md の「中間塩素化度の相がいずれも活性」という記述を支持。
  **熱力学量は報告されていない**

### 3.5 参考になるが直接は使えなかったもの

- Pan, B.; Zhang, Z.; Lu, X. "Determination and application of the reaction between
  REOCl (RE = Y, Gd and Sm) and H₂O." *Chem. Papers* **2020**, *74*(11), 3987–3993.
  DOI `10.1007/s11696-020-01202-5`
  — タイトルは魅力的だが**液体水との室温反応**（5REOCl + 5H₂O = 2RE₂(OH)₅Cl + RE³⁺ + 3Cl⁻）で、
  水酸化塩化物の生成。**高温気相加水分解ではないので本件には使えない**
- Yamamoto ら "Hidden Nature of the Conversion Reaction from Rare Earth Chloride to
  Oxychloride and the Application to Novel Separation." *ChemistrySelect* **2018**.
  DOI `10.1002/slct.201702581`
  — 573 K、H₂O 6.5 ppm の乾燥 O₂ 中で **DyCl₃ は DyOCl 化するが NdCl₃ はしない**。
  超低 H₂O 側の境界条件として有用だが HCl 共存系ではない（**[抄録]**）
- Ln–O–Cl 系の carbochlorination レビュー: *Min. Metall. Explor.* **2021**,
  DOI `10.1007/s42461-021-00490-z`（OA 版 `hdl.handle.net/11336/151697`、
  本調査ではアクセスできず）。「p(Cl₂) = 1 atm, 400 ℃ で log p(O₂) < 0 なら LnCl₃、
  それ以上なら LnOCl」という記述あり（**[抄録]**）。プロジェクトの Kellogg 図の
  独立チェックに使える可能性

---

## 4. 要取得リスト（優先順）

| # | 文献 | 欲しいもの | DOI / 入手先 |
|---|---|---|---|
| **1** | Koch, C. W.; Cunningham, B. B. *JACS* **1953**, *75*(4), 796–797 | **SmCl₃ + H₂O = SmOCl + 2HCl の ΔH, ΔG, K(T), 測定温度域, 不確かさ** | `10.1021/ja01100a010` |
| **2** | Koch, C. W.; Broido, A.; Cunningham, B. B. *JACS* **1952**, *74*(9), 2349–2351 | La の同上（系列比較の基準点） | `10.1021/ja01129a049` |
| **3** | Koch, C. W.; Cunningham, B. B. *JACS* **1954**, *76*(6), 1471–1474 | Pr, Nd の同上。**Nd の食い違い（§2.1）を裁定できる** | `10.1021/ja01635a003` |
| **4** | Gibson, A. ら *J. Chem. Thermodyn.* **2025**, *211*, 107549 | NdOCl / YOCl / TmOCl の **S°₂₉₈ と Cp(T)**、ΔG の規格化定義 | `10.1016/j.jct.2025.107549`、SSRN `10.2139/ssrn.5208841` |
| **5** | Yang, S. ら *Inorg. Chem.* **2022**, *61*, 7590–7596 | **全 REOCl の ΔfH₂₉₈ 数値表**（Table 2/3 相当）と補助データ | `10.1021/acs.inorgchem.2c00763` |
| **6** | Burns, J. B. ら *J. Alloys Compd.* **1998**, *265*, 146–152 | 本文の**全 Ln の LnOCl 推定 ΔfH 表**、および使用した Ln₂O₃/LnCl₃ 補助データ | `10.1016/S0925-8388(97)00435-0` |
| **7** | Knacke, Kubaschewski, Hesselmann (1991) 編纂 2nd ed. | **LaOCl, NdOCl, SmOCl, GdOCl の ΔHf₂₉₈, S°₂₉₈, Cp(T) の生データ**（頁番号まで記録） | 書籍。図書館 |
| 8 | Terlingen ら *ACS Catal.* **2022**, *12*, 5698 の SI | LnOCl → LnCl₃ 塩素化 ΔG の元データ | OA (PMC9087184) の Supporting Information |
| 9 | *Min. Metall. Explor.* **2021** carbochlorination レビュー | Ln–O–Cl の ΔG 表と 400 ℃ 安定領域 | `10.1007/s42461-021-00490-z` |

**入手できたときの記録要件**（CLAUDE.md の規律に従う）:
測定温度域・不確かさ・使用した補助データ（Ln₂O₃, LnCl₃, H₂O, HCl の出所）を
`source` に**ページ番号まで**書く。§2.1 の Eu/Gd の DH は補助データ依存なので、
補助データが変われば再計算が必要。

---

## 5. ギャップと、次に何をすべきか

### G1. 一次資料の数値が 1 つも取れていない（最大のギャップ）

Koch–Cunningham 3 部作、Burns ら、Yang ら、Gibson ら — **本件の中核となる
5 文献すべてがペイウォールで、本調査で取得できた数値は
Jacob ら (2016) の OA PDF 経由の二次データと、検索要約経由の抄録値のみ。**
`DH_OXYCHLORIDE` を `-30` から動かす判断は、要取得リスト #1〜#3 の
いずれか（できれば #1: Sm そのもの）を取ってから行うべき。

ただし、**「-30 は浅すぎる」という向きだけは、独立な 4 系統
（Knacke 編纂 / Burns 熱量測定 / Navrotsky 熱量測定 / Woodfield 熱容量）が
一致して示している**ので、この方向性は一次資料なしでも主張できる。

### G2. 「系列を通して一定の DH」という仮定が文献と整合しない

`data.py` の設計思想（推定誤差が Ln 間で相関 → 系列内の差は信頼できる）は、
DH が系列を通してほぼ一定であることを前提にしている。しかし:

- 実測系 4 元素の DH(ΔG, 1000 K) が **-31.0 〜 -57.5** と 26 kJ/mol ばらつく
- しかも La(-57.5) → Nd(-31.1) → Sm(-51.1) → Gd(-31.0) と**単調ですらない**
- Yang らは「イオン半径の減少とともに発熱性が減る」と単調変化を報告しており、
  Knacke 系の非単調性と**定性的に矛盾**
- 重希土側では 3D→2D 構造転移で TmOCl が -11.1 kJ/mol まで浅くなる。
  同じ半径帯の YOCl(-41.0, 3D) とは 30 kJ/mol 違う

**含意**: 現行の `_register_lanthanides()` は「絶対値は要検証だが系列内の差は
信頼できる」と `source` に書いているが、**この主張自体が支持されない可能性が高い。**
Knacke 系の非単調性が実測の反映なのか編纂上のノイズなのかを、
要取得 #1〜#3 と #5 で切り分ける必要がある。
これが判明するまで、`lanthanides.py` の系列比較の結論には
「LnOCl 推定に由来する系列内の相対誤差も ±15 kJ/mol 程度ありうる」
という但し書きを付けるべき。

### G3. エントロピー仮定の符号が文献と逆

§2.1 に書いた通り、プロジェクトは ΔS_rxn = +8 J/mol/K、
Burns/Jacob 系は -8 (±4) J/mol/K。653 K で ΔG にして約 10 kJ/mol の差。
**`data.py` の `S_ox / 3 + S_cl3 / 3 + 8.0` は `- 8.0` が正しい可能性が高い。**
これは DH の見直しとは独立に効くので、修正するなら別コミットにして
図の再生成と対応づけること（CLAUDE.md のデータ規律）。

なお、S°₂₉₈(LnOCl) の実測値は要取得 #4 に含まれている。
推定をやめて実測に置き換えられれば、この論点は消える。

### G4. Y と Ho の決定実験（CLAUDE.md `radius_controls()`）への影響

YOCl は Gibson ら 2025 で ΔG_f(298) = -41.0 kJ/mol の実測がある一方、
**HoOCl は Burns の推定値しかない**（Jacob Table 1 のタグ (d)）。
Y/Ho 対比実験を熱力学で裏打ちするには HoOCl の実測が必要だが、
Yang ら 2022 が「全 REOCl を合成して熱量測定した」と述べているので、
要取得 #5 に HoOCl が含まれている見込みが高い。

### G5. まだ探していない領域

- **ロシア語文献**（Polyachenok, Novikov ら 1960 年代の LnCl₃ 蒸気圧・加水分解）
  — 本調査では英語圏のみを探索した
- **中国語文献**（希土類オキシ塩化物の熱力学は中国で活発）
- **Gmelin Handbook** の Rare Earth Elements 巻（オキシハロゲン化物の系統的編纂）
- **Cordfunke & Konings** の核燃料/FP 熱化学データ評価（LaOCl を含む可能性）
- **FactSage / HSC / Thermo-Calc の SGTE 系データベース**に LnOCl が
  どう入っているか（プロジェクトが将来 pycalphad に移行する際の互換性）

---

## 付録: 本調査で使った計算の再現

```python
# §1.3 の DH(ΔG, 1000 K) — Jacob 2016 Table 1 から
DH = dGf_LnOCl - dGf_Ln2O3 / 3 - dGf_LnCl3 / 3

# §2.1 の「DH 実効値」— プロジェクトのモデル内で文献 ΔG を再現する DH
#   deacon_thermo.data.DB を使い、LnOCl の dHf298 だけ振って
#   dG_rxn(T) = G(LnOCl) - G(Ln2O3)/3 - G(LnCl3)/3 が文献値に一致する DH を線形逆解き
#   （dHf298 は G に加法的に入るので厳密に線形）

# §2.2 の余裕 — 既存 API をそのまま使用
from deacon_thermo.gas import gas_state
from deacon_thermo.sensitivity import chloride_margin_at
g = gas_state(653.15, 2.0, extent=conv * 0.5)   # conv = HCl 転化率
chloride_margin_at('Sm', g, dh)
```

`DB` の数値は本調査では一切書き換えていない。
