# Cu 塩化物気相種の熱力学データ調査

調査日: 2026-08-05
対象: `src/deacon_thermo/data.py` の CuCl(g), Cu3Cl3(g), CuCl2(g) エントリの一次資料照合
結論を先に書く: **Cu3Cl3(g) の暫定値は一次資料（NIST-JANAF）と大きく食い違っており、揮発計算に致命的な影響がある。CuCl(g) は一次資料と一致し妥当。CuCl2(g) は一次資料が見つからず未決着。**

---

## 1. サマリ表

| 種 | 現行暫定値 ΔHf₂₉₈ [kJ/mol] | 一次資料値 | 現行暫定値 S°₂₉₈ [J/mol/K] | 一次資料値 | 信頼度 | 出典 |
|---|---|---|---|---|---|---|
| CuCl(g) | 91.1 | **91.086** | 237.0 | **237.207** | **全文確認（一致）** | NIST-JANAF 4th ed. (Chase 1998), table Cl-012、review date 1966 |
| Cu3Cl3(g) | -305.0 | **-258.571** | 469.0 | **429.526** | **全文確認（不一致・大）** | NIST-JANAF 4th ed. (Chase 1998), table Cl-132、review date 1966 |
| CuCl2(g) | -43.0 | **見つからず** | 278.0 | **見つからず** | 未確認 | — |

(参考: JANAF は Cp₂₉₈ も与える。Cu3Cl3(g): 124.570 J/mol/K, CuCl(g): 35.261 J/mol/K)

NIST-JANAF は NIST WebBook 経由（janaf.nist.gov のテーブル HTML/TXT を直接 2 回独立に取得し、数値が一致することを確認済み）と NIST Chemistry WebBook（webbook.nist.gov、同じ Chase 1998 を出典として引用）の両方から同一の値が得られており、**CuCl(g) と Cu3Cl3(g) については "全文確認" と言ってよい信頼度**。

副産物として、計算に絡む凝縮相も確認できた:

| 種 | 現行暫定値 ΔHf₂₉₈ | JANAF 値 | 現行暫定値 S°₂₉₈ | JANAF 値 |
|---|---|---|---|---|
| CuCl(s) | -137.2 (Barin) | -138.072 | 86.2 | 87.027 |
| CuCl2(s) | -220.1 (Barin) | **-205.853** | 108.1 | 108.085 |

CuCl(s) はほぼ一致（誤差 1 kJ/mol 未満）。**CuCl2(s) は ΔHf が 14.2 kJ/mol も食い違っている**（エントロピーはほぼ一致）。これは今回の依頼スコープ外だが、`volatility.py` の `partial_pressures()` が CuCl2(g) の基準相として CuCl2(s) を使っているため、CuCl2(g) の評価に間接的に影響する。詳細は §5 のギャップ参照。

---

## 2. Cu3Cl3(g): 現行値との差分、および 653 K での蒸気圧への影響

`volatility.py` の `partial_pressures()` は次を計算している（コード確認済み）:

```
3 CuCl(s) -> Cu3Cl3(g)
dG(T) = G(Cu3Cl3,g,T) - 3 * G(CuCl, supercooled liquid, T)
p(Cu3Cl3) ∝ exp(-dG / RT)
```

この反応の ΔH, ΔS を現行暫定値と JANAF 値で比較する（CuCl(s) は各データセット内で自己無矛盾になるよう対応する値を使用）。

**現行暫定値（Cu3Cl3(g): -305.0 / 469.0, CuCl(s) 現行値: -137.2 / 86.2 [Barin]）**

- ΔH_rxn = -305.0 - 3×(-137.2) = **+106.6 kJ/mol**
- ΔS_rxn = 469.0 - 3×86.2 = **+210.4 J/mol/K**
- ΔG_rxn(653 K) = 106.6 - 653×0.2104 = **-30.8 kJ/mol**（負 = 反応は自発的に進む）

**JANAF 一次資料値（Cu3Cl3(g): -258.571 / 429.526, CuCl(s): -138.072 / 87.027）**

- ΔH_rxn = -258.571 - 3×(-138.072) = **+155.6 kJ/mol**
- ΔS_rxn = 429.526 - 3×87.027 = **+168.4 J/mol/K**
- ΔG_rxn(653 K) = 155.6 - 653×0.1684 = **+45.6 kJ/mol**（正 = 反応は不利）

**差**: ΔΔG = 45.6 - (-30.8) = **+76.4 kJ/mol**（653 K）

蒸気圧比 = exp(-ΔΔG / RT)、RT(653 K) = 8.314×653/1000 = 5.429 kJ/mol

```
exp(-76.4 / 5.429) = exp(-14.07) ≈ 7.7×10⁻⁷
```

**つまり、JANAF の一次資料値を使うと、現行暫定値を使った場合に比べて p(Cu3Cl3) の予測値がおよそ 10⁻⁶〜10⁻⁷ 倍（約 130 万分の 1）になる。** 現行の -305.0/469.0 という値は Cu3Cl3(g) を実際より大幅に安定・高エントロピーに見積もっており、その結果 Cu 揮発速度を桁違いに過大評価している可能性が高い。

これは `CLAUDE.md` に記載されている「理想 Temkin モデルは実測寿命を約 2 桁外す（寿命を過小評価する＝揮発を過大評価する）」という既知の問題と**同じ方向**の効果であり、Cu3Cl3(g) データの誤りだけで 2 桁どころか 6 桁分の変化を生みうる規模になる。クロロ銅酸錯体による安定化を考慮する前に、まずこの気相種データそのものを JANAF 値に差し替えて再計算する価値が高い。

**参考（未検証の傍証、要一次資料確認）**: Guido, Balducci, Gigli, Spoliti, *J. Chem. Phys.* **55**, 4566 (1971) の Knudsen 効果流出質量分析法による研究が、3CuCl(s) → Cu3Cl3(g) の蒸発エンタルピー ΔvapH°(640 K) ≈ 141 kJ/mol を報告しているという言及を検索エンジンの要約経由で見つけた。この値は JANAF から導かれる ΔH_rxn(298K)=155.6 kJ/mol と同程度のオーダーで整合的（温度差・Cp補正を考えれば矛盾しない）であり、JANAF 値を支持する方向の傍証になる。**ただし原論文本文を直接確認できていないため、この数値自体は「検索要約経由・未確認」として扱うこと。** 要取得リストに追加した。

---

## 3. CuCl(g): 現行値の妥当性

現行値 ΔHf₂₉₈ = 91.1 kJ/mol, S°₂₉₈ = 237.0 J/mol/K は NIST-JANAF 4th ed.（Chase 1998, review date 1966）の値 91.086 / 237.207 とほぼ完全に一致した（小数点以下の丸めのみの差）。data.py の出典欄は "Barin" だが、Barin の値も同じ 1966 年の一次データ系譜（JANAF/NBS の評価）に由来しているとみられ、実質的に同じ数値である。

**結論: CuCl(g) の現行値は妥当。信頼度を FAIR → GOOD に上げてよい（出典を JANAF に更新することを推奨、`data.py` の変更は今回のスコープ外なので提案のみ）。**

---

## 4. CuCl2(g): 見つからなかったこと

- **NIST-JANAF**: `Cu-index.html` を直接取得して確認したが、Cu-Cl 種として掲載されているのは CuCl(cr, l, cr-l, g) と Cu3Cl3(g) のみ。**CuCl2 は結晶相（Cl-079.html）しか JANAF に存在しない。CuCl2(g) の JANAF テーブルは存在しない。**
- **NIST Chemistry WebBook**（webbook.nist.gov, ID=C7447394、copper dichloride）: 固相の熱化学データ（JANAF 由来、上記 §1 参照）とイオンエネルギー論データ（イオン化エネルギー等）はあるが、**気相の ΔfH°・S° は掲載されていない。**
- **Barin, Thermochemical Data of Pure Substances**: 索引には "CuCl2" が掲載されている（p.610 付近）が、同じ索引で CuCl[g] と Cu3Cl3[g] には明示的に "(GAS)" と注記されているのに対し、CuCl2 にはその注記がない。これは Barin でも CuCl2 が凝縮相のみ収録で、気相データがない可能性を示唆する（ただし実データページ本文までは今回アクセスできておらず、確定はできない）。
- **分光学的存在確認**: CuCl2(g) という分子種自体は実在が確認されている。
  - Hougen & Leroi, *J. Chem. Phys.* **34**, 1670 (1961), "Application of Ligand Field Theory to the Electronic Spectra of Gaseous CuCl2, NiCl2, and CoCl2" — 1000 ℃ 前後で気相 CuCl2 の電子吸収スペクトルを測定。D∞h 対称の直線分子として議論。
  - (著者未特定) *J. Chem. Phys.* **44**, 4387 (1966), "Electronic Absorption Spectra of the Gaseous 3d Transition-Metal Dichlorides"(VCl2, CrCl2, FeCl2, CoCl2, NiCl2, CuCl2 を含む) — 同様に気相存在を確認。
  - これらは分子構造・振動スペクトルの情報であり、**ΔHf・S° の値そのものは含まれていない**（statistical mechanics による third-law 計算をすれば S° は導出できる可能性があるが、それには振動振動数・回転定数・電子基底状態の縮重度などが必要で、今回の調査range外）。
- **Schäfer 系のガス輸送反応研究**: Wächter & Schäfer, *Z. Anorg. Allg. Chem.* (1980), "Das Gleichgewicht CuCl2,f + Al2Cl6,g = CuAl2Cl8,g und die beteiligten Nebenreaktionen" というタイトルの論文の存在を確認した。CuCl2 固相と Al2Cl6(g) の輸送平衡を扱っており、"begleitende Nebenreaktionen"（随伴する副反応）に CuCl2(g) や CuCl(g)+1/2Cl2(g) 系の平衡が含まれている可能性がある。**本文未確認、要取得。**

**結論: CuCl2(g) の ΔHf₂₉₈・S°₂₉₈ は一次資料で確認も否定もできなかった。現行の -43.0 / 278.0 kJ, J/mol/K は "要検証" のまま。** data.py の `EST` (ESTIMATE) 信頼度は正しい格付けであり、変更すべきでない。

---

## 5. 出典詳細

### 全文確認（信頼度: 高）

1. **NIST-JANAF Thermochemical Tables, 4th ed.** (M.W. Chase Jr., 1998, NIST Standard Reference Database 13)
   - CuCl(g): table `Cl-012` (janaf.nist.gov/tables/Cl-012.txt) — 2 回独立に取得・値一致確認済み
   - Cu3Cl3(g)（"Cl3Cu3"表記）: table `Cl-132` (janaf.nist.gov/tables/Cl-132.html および .txt) — 2 回独立に取得・値一致確認済み
   - CuCl(cr): table `Cl-009` (janaf.nist.gov/tables/Cl-009.txt)
   - CuCl2(cr): table `Cl-079` (janaf.nist.gov/tables/Cl-079.txt)
   - 索引: janaf.nist.gov/tables/Cu-index.html（Cu-Cl 種の一覧を確認。CuCl2(g), Cu4Cl4(g) 等は掲載なし）
   - いずれも "Enthalpy Reference Temperature Tr = 298.15 K", "Standard State Pressure p° = 0.1 MPa", "data reviewed March 1966" と明記されている。

2. **NIST Chemistry WebBook** (webbook.nist.gov)
   - copper chloride CuCl: ID=C7758896 — 上記 JANAF と同一値を掲載、出典として Chase 1998 を明記
   - copper dichloride CuCl2: ID=C7447394 — 固相データのみ、出典 Chase 1998
   - tricopper trichloride Cu3Cl3: ID=C11093655 — イオンエネルギー論データのみ掲載、熱化学（ΔHf, S°）テーブルは無し

### 二次資料経由・書誌確認のみ（本文未取得）

3. **Barin, I., *Thermochemical Data of Pure Substances*, 3rd ed. (VCH, 1995)** — 索引ページ（pdfcoffee.com 掲載の電子版目次）で CuCl[g], CuCl2, Cu3Cl3[g] がいずれも p.610 付近に掲載されていることを確認したが、数値本文は未取得。data.py の CuCl(g)/CuCl(s)/CuCl2(s) の出典表記 "Barin" の妥当性検証はできたが、CuCl2(g) 相当のデータが実在するかは未確定。

4. **Guido, M.; Balducci, G.; Gigli, G.; Spoliti, M.** "Mass spectrometric study of the vaporization of cuprous chloride and the dissociation energy of Cu3Cl3, Cu4Cl4, and Cu5Cl5", *J. Chem. Phys.* **55**, 4566 (1971). — Knudsen 効果流出質量分析。CuCl(s) 蒸気中に Cu3Cl3, Cu4Cl4 が同程度の量、Cu5Cl5 が微量存在すると報告（検索エンジン要約より）。ΔvapH°(640K) ≈ 141 kJ/mol という数値の言及があったが**原文未確認、要検証**。

5. **Guido, M.; Gigli, G.; Balducci, G.** "Dissociation Energy of CuCl and Cu2Cl2 Gaseous Molecules", *J. Chem. Phys.* **57**, 3731 (1972). — CuCl・Cu2Cl2 の解離エネルギー。本文未取得。

6. **Hougen, J.T.; Leroi, G.E.** "Application of Ligand Field Theory to the Electronic Spectra of Gaseous CuCl2, NiCl2, and CoCl2", *J. Chem. Phys.* **34**, 1670 (1961). — 気相 CuCl2 の分光学的実在確認。ΔHf/S° は含まない。本文未取得。

7. "Electronic Absorption Spectra of the Gaseous 3d Transition-Metal Dichlorides", *J. Chem. Phys.* **44**, 4387 (1966). — 著者名は検索結果から特定できず。CuCl2(g) を含む 3d 金属ジクロリド気相種の吸収スペクトル。本文未取得。

8. **Wächter, H.; Schäfer, H.** "Das Gleichgewicht CuCl2,f + Al2Cl6,g = CuAl2Cl8,g und die beteiligten Nebenreaktionen", *Z. Anorg. Allg. Chem.* (1980). — CuCl2 固相と Al2Cl6 気相の輸送平衡研究。CuCl2(g) 関連データを含む可能性があるが本文未取得。書誌の巻・ページ番号も未確定（要取得）。

9. **Schäfer, H.; Binnewies, M.** "Die Stabilität gasförmiger Dimerer Chloridmolekeln", *Z. Anorg. Allg. Chem.* (1974). — 気相二量体塩化物分子の安定性の一般論。Cu2Cl2(g) に触れている可能性あり。本文未取得。

---

## 6. 要取得リスト（優先順）

1. **Barin (1995) 実データページ本文（p.610 前後）** — CuCl2 の項目が本当に気相データを含まないか、含む場合はその数値を確認。全文書誌: Barin, I. *Thermochemical Data of Pure Substances*, 3rd ed., VCH Verlagsgesellschaft, Weinheim, 1995. ISBN 978-3-527-28531-0.
2. **Guido, Balducci, Gigli, Spoliti (1971) *J. Chem. Phys.* 55, 4566** — DOI: 10.1063/1.1676789（要確認）。Cu3Cl3(g) の JANAF 値の独立検証、および ΔvapH°(640K)≈141 kJ/mol という数値の真偽確認に必須。
3. **Wächter & Schäfer (1980) *Z. Anorg. Allg. Chem.*** — 巻・ページ番号を含む完全書誌の確定と、CuCl2(g) 関連データの有無確認。CuCl2(g) の唯一の手がかりになりうる。
4. **Guido, Gigli, Balducci (1972) *J. Chem. Phys.* 57, 3731** — DOI: 10.1063/1.1678845（要確認）。CuCl/Cu2Cl2 解離エネルギーとJANAF値の相互検証。
5. Gurvich, L.V. ら (IVTANTHERMO) の評価値 — 今回アクセスできなかった。ロシア語圏の評価が JANAF と独立クロスチェックになる可能性がある。

## 7. ギャップ（今後の方針への示唆）

- **最優先で `data.py` の Cu3Cl3(g) を JANAF 値（ΔHf₂₉₈ = -258.571 kJ/mol, S°₂₉₈ = 429.526 J/mol/K, source="JANAF Cl-132"）に差し替えることを推奨する。** §2 の計算通り、現行値との差は 653 K で蒸気圧にして 6〜7 桁に達し、Cu 揮発の結論を直接左右する。CLAUDE.md が「理想 Temkin モデルが実測寿命を 2 桁外す」としている問題の一部（あるいは大部分）が、クロロ銅酸錯体の安定化ではなく単純に Cu3Cl3(g) の暫定値の誤りで説明できる可能性がある。差し替え後に `test_ideal_model_reproduces_observed_lifetime` の xfail 具合がどう変わるか確認する価値が高い。
- CuCl(g) は現行値のまま出典を "Barin" → "JANAF (Cl-012, Chase 1998)" に更新し、信頼度を FAIR → GOOD にしてよい。
- CuCl2(g) は今回の調査では一次資料が見つからなかった。EST のまま据え置くのが適切。Wächter & Schäfer (1980) の入手が次の一手。
- 副産物として見つかった CuCl2(s) の ΔHf 不一致（現行 -220.1 kJ/mol vs JANAF -205.853 kJ/mol、差 14.2 kJ/mol）は、`volatility.py` が CuCl2(g) の基準相として CuCl2(s) を使う設計上、CuCl2(g) 評価にも波及する。CuCl2(g) データが得られた際は CuCl2(s) も同時に JANAF 値へ差し替えて整合性を取ること。
