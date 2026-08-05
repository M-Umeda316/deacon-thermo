# アルカリ–ランタノイド／アルカリ–銅 塩化物複塩の熱力学データ 文献調査

調査日: 2026-08-05
対象: Cu-K-Ln/γ-Al₂O₃ Deacon 触媒（380 ℃ = 653 K, HCl/O₂ = 2:1）の凝縮相安定性・
Cu 活量低下効果の評価に必要な複塩生成熱力学

> **信頼度ラベルの定義**
> - **【全文確認】** … 論文全文（PDF/HTML）を実際に取得して読んだ数値。表番号・式番号を併記。
> - **【二次資料経由】** … 別の論文・データベースが引用している値、または CALPHAD 最適化により
>   導かれたモデルパラメータ。原著の一次データではない。
> - **【抄録のみ】** … 抄録・要旨・検索スニペットのみから取った記述。数値の桁・符号は
>   原典で必ず再確認すること。
>
> **本レポートに書かれていない数値は「見つからなかった」ということ。** 推定値・記憶による
> 値は一切書いていない。

---

## 0. 結論の要約（先に読む用）

1. **K–Ln 複塩（K₃LnCl₆ / K₂LnCl₅ / KLn₂Cl₇）は 380 ℃ ではすべて固体。**
   KCl–LnCl₃ 二元系の最低共晶・最低不変反応温度は Ln = Nd で 719 K（446 ℃）、
   Pr で 762 K（489 ℃）であり、653 K を上回る。純粋な K–Ln 塩化物系だけでは
   反応温度で液相は現れない。
2. **一方 K–Cu 塩化物は 380 ℃ で完全に液体。** KCl–CuCl 共晶は 145.9 ℃、
   KCl–CuCl₂ 系は 533–543 K（260–270 ℃）以降に相変化が連続する。
   つまり **実際の触媒相は「Cu 塩化物融液に KCl と LnCl₃ が溶けたもの」** であって、
   K–Ln 複塩が独立の固相として析出するかどうかは Cu 融液中の活量次第。
   → `melt.py` の活量モデルでは K–Ln 複塩を「固相候補」として扱うのが妥当。
3. **K₃LnCl₆ の低温安定限界がランタノイド収縮とともに劇的に下がる**（本レポートの
   最重要な系列トレンド）。La では K₃LaCl₆ 自体が存在せず、Pr は 489 ℃ 以上、
   Nd は 446 ℃ 以上、Dy では **39 ℃** 以上で安定。Sm はこの間にある（要一次資料）。
4. **CuCl–CuCl₂–KCl 三元融液の平衡塩素圧データが実在する**
   （Fontana et al., Ind. Eng. Chem. 44 (1952) 363）。これは CLAUDE.md の
   作業優先順位 2「`fit_interactions()` の較正」に直接使える一次資料。**最優先で入手すべき。**
5. **数値として確定できたのは Pr 系 2 化合物・Nd 系 3 化合物の計 5 件のみ**
   （いずれも CALPHAD 最適化値＝二次資料）。**溶解熱量測定による一次データ（Seifert 系列、
   Blachnik & Selle）はすべてペイウォールで取得できなかった。**

---

## 1. サマリ表

### 1-1. K–Ln 複塩：二元塩化物からの生成量

反応定義（すべて固体、298.15 K 基準の標準状態を採るが、CALPHAD 値は
Neumann–Kopp 則（ΔC_p = 0）を仮定しているため全温度域で同一の ΔH, ΔS を与える）:

- (A) `3 KCl(s) + LnCl₃(s) → K₃LnCl₆(s)`
- (B) `2 KCl(s) + LnCl₃(s) → K₂LnCl₅(s)`
- (C) `KCl(s) + 2 LnCl₃(s) → KLn₂Cl₇(s)`
- (D) シンプロポーショネーション（Seifert の "Synreaktion"）: `KCl(s) + K₂LnCl₅(s) → K₃LnCl₆(s)`

| 化合物 | 反応 | ΔH / kJ mol⁻¹ | ΔS / J mol⁻¹ K⁻¹ | 融点・分解温度 | 380 ℃での状態 | 信頼度 | 出典 |
|---|---|---|---|---|---|---|---|
| K₂LaCl₅ | (B) | **未取得** | 未取得 | 913 K 一致融解（Seifert）／916 K（Song） | 固体 | 【二次資料経由】温度のみ | Gong et al. arXiv:2406.15223 Table 5、原典 Seifert 1985 [S-La], Blachnik 1971 [B-La] |
| K₃La₅Cl₁₈ | — | **未取得** | 未取得 | 885 K 包晶 `L + LaCl₃ → K₃La₅Cl₁₈` | 固体 | 【二次資料経由】 | 同上 |
| K₃LaCl₆ | (A) | — | — | **存在しない**（KCl–LaCl₃ 系に生成相なし） | — | 【二次資料経由】 | 同上 |
| K₂PrCl₅ | (B) | **−50.387** | 0（温度非依存項なし） | 890 K 包晶生成、室温まで安定 | 固体 | 【全文確認】ただし CALPHAD 最適化値 | Sridar, Hao, Xiong, Calphad 81 (2023) 102552, Table 2 |
| K₃PrCl₆ | (A) | **+11.199** | **+80.774** | **762 K 以上でのみ存在**、945 K 一致融解 | 固体（安定域外＝KCl+K₂PrCl₅ に分解） | 【全文確認】CALPHAD 最適化値 | 同上, Table 2, Table 3 |
| K₃PrCl₆ | (D)（上2行から導出） | +61.586 | +80.774 | 転移温度 762.5 K | — | 【二次資料経由】導出値 | 同上 |
| K₂NdCl₅ | (B) | **−32.900** | **+18.57** | 883 K 包晶生成、室温まで安定 | 固体 | 【全文確認】CALPHAD 最適化値 | Hao, Sridar, Xiong, J. Mol. Liq./OSTI 2259217, Table 1 |
| K₃NdCl₆ | (A) | **+17.097** | **+88.5** | **719 K 以上でのみ存在**、972 K 一致融解 | 固体（安定域外） | 【全文確認】CALPHAD 最適化値 | 同上 |
| KNd₂Cl₇ | (C) | **+6.313** | **+39.91** | 742 K 以上でのみ存在、783 K 包晶分解 | 固体（安定域外） | 【全文確認】CALPHAD 最適化値 | 同上 |
| K₃NdCl₆ | (D)（上2行から導出） | +49.997 | +69.93 | 転移温度 715 K（実験値 719 K） | — | 【二次資料経由】導出値 | 同上 |
| K₃LnCl₆ (Ln = Ce, Pr, Nd) | (D) | **+48 〜 +55** | **+63 〜 +65** | いずれも高温相 | — | 【抄録のみ】 | Seifert, J. Therm. Anal. Cal. 67 (2002) 789 / 83 (2006) 479 の抄録；Gaune-Escard & Rycerz |
| K₃TbCl₆ | — | — | 固–固転移 ΔtrsH = **6.1 kJ mol⁻¹** | 一致融解（温度は本文中、未取得） | 固体 | 【抄録のみ】 | Rycerz & Gaune-Escard, JTAC 68 (2002) 973 抄録 |
| KTb₂Cl₇ | — | — | K 塩の固–固転移値は抄録に記載なし（Rb: 17.1, Cs: 12.1/10.9 kJ mol⁻¹） | 一致融解 | 固体 | 【抄録のみ】 | 同上 |
| K₃DyCl₆ | (D) | 正（格子エンタルピーの損失）／数値未取得 | 正／数値未取得 | **312 K（39 ℃）以上でエントロピー安定化** | 固体（安定） | 【抄録のみ】 | Seifert & Krämer, ZAAC 620 (1994) 1543 抄録 |
| K₂DyCl₅ | (B) | 数値未取得（溶解熱量測定済み） | — | K₂PrCl₅ 型構造 | 固体 | 【抄録のみ】 | 同上 |
| KDy₂Cl₇ | (C) | 数値未取得（溶解熱量測定済み） | — | — | 固体 | 【抄録のみ】 | 同上 |
| K₂EuCl₅ | (B) | 数値未取得（溶解熱量＋EMF 測定済み） | — | KCl 系で室温安定なのは A₂EuCl₅ のみ | 固体 | 【抄録のみ】 | Seifert & Sandrock, ZAAC 587 (1990) 110 抄録 |
| K₃GdCl₆ / K₂GdCl₅ / KGd₂Cl₇ | (A)(B)(C) | 数値未取得（溶解熱量＋EMF 測定済み） | — | K₂GdCl₅ は K₂PrCl₅ 型 | 固体 | 【抄録のみ】 | Seifert, Sandrock & Thiel, ZAAC 598 (1991) 307 抄録 |
| **K₃SmCl₆ / K₂SmCl₅ / KSm₂Cl₇** | (A)(B)(C) | **数値未取得（本プロジェクトの参照系そのもの）** | — | — | 固体（推定） | **未入手** | Thiel & Seifert, Thermochim. Acta 133 (1988) 275 |

**KCl–LnCl₃ 二元系の不変反応（実験値）**

| 系 | 反応 | T / K | x(LnCl₃) | 信頼度 | 出典 |
|---|---|---|---|---|---|
| KCl–LaCl₃ | `L → KCl + K₂LaCl₅` | 853（Seifert）／845（Song） | 0.22 | 【二次資料経由】 | arXiv:2406.15223 Table 5 |
| KCl–LaCl₃ | `L → K₂LaCl₅`（一致融解） | 913／916 | 0.333 | 同上 | 同上 |
| KCl–LaCl₃ | `L → K₂LaCl₅ + K₃La₅Cl₁₈` | 851 | 0.51 | 同上 | 同上 |
| KCl–LaCl₃ | `L + LaCl₃ → K₃La₅Cl₁₈` | 885 | — | 同上 | 同上 |
| KCl–PrCl₃ | `L → KCl + K₃PrCl₆` | 886 | 0.165 | 【全文確認】(実験値は Seifert 1987) | Calphad 81 (2023) 102552 Table 3 |
| KCl–PrCl₃ | `L → K₃PrCl₆`（一致融解） | 945（944） | 0.25 | 同上 | 同上 |
| KCl–PrCl₃ | `L + K₃PrCl₆ → K₂PrCl₅` | 890 | 0.345 | 同上 | 同上 |
| KCl–PrCl₃ | `K₃PrCl₆ → KCl + K₂PrCl₅` | 762（768） | — | 同上 | 同上 |
| KCl–PrCl₃ | `L → PrCl₃ + K₂PrCl₅` | 772 | 0.560 | 同上 | 同上 |
| KCl–NdCl₃ | `L → KCl + K₃NdCl₆` | 899 | 0.160 | 【全文確認】(実験値は Seifert 1988) | OSTI 2259217 Table 4 |
| KCl–NdCl₃ | `L → K₃NdCl₆`（一致融解） | 972 | 0.25 | 同上 | 同上 |
| KCl–NdCl₃ | `L + K₃NdCl₆ → K₂NdCl₅` | 883 | 0.335 | 同上 | 同上 |
| KCl–NdCl₃ | `L + NdCl₃ → KNd₂Cl₇` | 783 | 0.600 | 同上 | 同上 |
| KCl–NdCl₃ | `L → K₂NdCl₅ + KNd₂Cl₇` | 766 | 0.540 | 同上 | 同上 |
| KCl–NdCl₃ | `KNd₂Cl₇ → K₂NdCl₅ + NdCl₃` | 742 | — | 同上 | 同上 |
| KCl–NdCl₃ | `K₃NdCl₆ → KCl + K₂NdCl₅` | 719 | — | 同上 | 同上 |

**液相の会合種（Cu 融液モデルへの参考）**

- `3 KCl(l) + NdCl₃(l) → K₃NdCl₆(l)`（液相中の中性会合種）: **ΔG = −73 900 J mol⁻¹**（温度非依存）
  【全文確認・CALPHAD 最適化値】OSTI 2259217 Table 1。
  → KCl–LnCl₃ 融液の短距離秩序（SRO）は `K₃LnCl₆` 会合体で表現できる、という
  モデル上の事実。Cu 融液に K/Ln を入れたときの活量低下を Temkin＋会合種で
  書くときの出発点になる。
- KCl–LaCl₃ 液相の混合エンタルピー実測: x(LaCl₃) = 0.496 で **−15.319 kJ mol⁻¹**
  （1173 K, Papatheodorou & Ostvold 1974）【二次資料経由】arXiv:2406.15223。

### 1-2. K–Cu 複塩

反応定義: (E) `KCl(s) + CuCl₂(s) → KCuCl₃(s)`、(F) `2 KCl(s) + CuCl(s) → K₂CuCl₃(s)`

| 化合物・不変点 | 系 | 値 | 380 ℃での状態 | 信頼度 | 出典 |
|---|---|---|---|---|---|
| K₂CuCl₃ | KCl–CuCl | 化学量論化合物として扱われる。ΔH_r は ab initio 計算値を固定、S° を最適化。**数値は未取得** | 液体 | 【抄録のみ】 | Niazi et al., Materialia 21 (2022) 101296 |
| KCl–CuCl 共晶 | KCl–CuCl | **145.9 ℃（419.0 K）, 64.9 mol% CuCl**（CALPHAD 計算値） | 液体 | 【抄録のみ】 | 同上 |
| KCl–CuCl 包晶 | KCl–CuCl | **241.2 ℃（514.4 K）** | 液体 | 【抄録のみ】 | 同上 |
| KCl–CuCl 熱安定性 | — | 500 ℃ まで質量減少なし。加水分解試験で HCl 検出されず | — | 【抄録のみ】 | 同上 |
| KCuCl₃ | KCl–CuCl₂ | 無水 KCl/CuCl₂ 状態図では**一致融解化合物**。融点の数値は未取得 | 液体（融液中） | 【二次資料経由】Wikipedia 経由 | — |
| K₂CuCl₄（無水） | KCl–CuCl₂ | **無水物は状態図に現れない。**加熱すると KCl·CuCl₂ + KCuCl₃ の混合物に変わる | — | 【二次資料経由】Wikipedia 経由 | — |
| K₂CuCl₄·2H₂O | — | **93 ℃ 以上で KCl + KCuCl₃ + H₂O に分解** | — | 【二次資料経由】Wikipedia 経由 | — |
| CuCl₂–KCl の相変化 | KCl–CuCl₂ | DTA で **533–543 K, 573–593 K, 603–623 K** に相変化（再現性を議論） | 液体を含む | 【抄録のみ】 | Żurowski, J. Therm. Anal. 45 (1995) 437 |
| CuCl₂–KCl の構成 | KCl–CuCl₂ | Cu/K モル比 0.2–2.0 を空気/Ar 下で TG。相変化と系の構成を議論 | — | 【抄録のみ】 | Żurowski, J. Therm. Anal. 44 (1995) 197 |
| KCuCl₃（DFT） | — | 生成エネルギー **−1.318 eV/atom**（元素基準、OQMD の基底状態） | — | 【検索スニペットのみ・未検証】 | oqmd.org/materials/composition/KCuCl3（調査時 502 エラーで再確認できず） |

> **KCuCl₃ / K₂CuCl₄ の生成エンタルピー（二元塩化物からの ΔH）は、
> 実験値・評価値ともに見つからなかった。** これは本調査の最大のギャップ。

---

## 2. Ln 系列トレンドの考察

本プロジェクトの方針（「系列内の差は絶対値より信頼できる」）に沿って、
絶対値ではなく **La → Er でどちらに動くか** を整理する。

### 2-1. トレンド A：K₃LnCl₆ の低温安定限界が急激に下がる（最も明瞭）

`K₃LnCl₆ ⇌ KCl + K₂LnCl₅` の平衡温度（これ以下では K₃LnCl₆ が分解）:

| Ln | Ln³⁺ 半径 (CN=6) / Å | K₃LnCl₆ 分解温度 | 出典・信頼度 |
|---|---|---|---|
| La | 1.032 | **化合物が存在しない** | 【二次資料経由】arXiv:2406.15223 |
| Ce | 1.01 | 高温相のみ（温度値未取得） | 【抄録のみ】Seifert 1986 |
| Pr | 0.99 | **762 K（489 ℃）** | 【全文確認】Calphad 81 (2023) |
| Nd | 0.983 | **719 K（446 ℃）** | 【全文確認】OSTI 2259217 |
| Sm | 0.958 | **未取得（要 Thiel & Seifert 1988）** | — |
| Eu | 0.947 | 高温相のみ（K 系で室温安定なのは K₂EuCl₅） | 【抄録のみ】ZAAC 587 (1990) |
| Gd | 0.938 | Sm/Eu 系と「対応する」と記載、数値未取得 | 【抄録のみ】ZAAC 598 (1991) |
| Dy | 0.912 | **312 K（39 ℃）** | 【抄録のみ】ZAAC 620 (1994) |
| Ho, Er, Y | 0.901 / 0.890 / 0.900 | 全系で K₃LnCl₆ が存在（Ce–Lu, Y）、温度未取得 | 【抄録のみ】ZAAC 627 (2001) |

> **半径 1.032 → 0.912 Å（約 12 %）の収縮で、分解温度が「存在しない」から
> 489 ℃ → 446 ℃ → 39 ℃ へ落ちる。** 380 ℃（653 K）を跨ぐのは
> **Nd（446 ℃）と Dy（39 ℃）の間**、すなわち **Sm–Gd 付近**。
> 参照系の Sm がまさにこの境界にある。SmCl₃ 25–30 wt% の系で、
> K₃SmCl₆ が固相として存在しうるか否かは 380 ℃ でぎりぎりの判定になる。
> → **Thiel & Seifert, Thermochim. Acta 133 (1988) 275–282 の入手が最優先。**

なお、K₃LnCl₆ が高温相である熱力学的理由は Seifert が明確に述べている:
`KCl + K₂LnCl₅ → K₃LnCl₆` は **格子エネルギーを失う（ΔH > 0）反応** であり、
K₃LnCl₆ が孤立配位多面体（[LnCl₆]³⁻）をもつためエントロピーが大きく増える。
TΔS > ΔH となる温度以上でのみ安定化する。Ln が小さくなると [LnCl₆]³⁻ が
安定になり ΔH が小さくなるため、転移温度が下がる、という描像で一貫している。

### 2-2. トレンド B：K₂LnCl₅ は系列全体で室温安定、しかし ΔH は Ln で大きく変わる

- Pr: `2KCl + PrCl₃ → K₂PrCl₅` ΔH = −50.4 kJ/mol
- Nd: `2KCl + NdCl₃ → K₂NdCl₅` ΔH = −32.9 kJ/mol

**この差（17.5 kJ/mol）を系列トレンドとして読んではいけない。**
両者は別々の CALPHAD 評価から来ており、液相モデルが異なる（Pr: ionic two-sublattice
＋ Redlich-Kister、Nd: ionic two-sublattice ＋ 中性会合種 K₃NdCl₆）。さらに
Pr の評価では **Seifert の実測生成エンタルピーを「相図が再現できない」という理由で
意図的に除外している**（Calphad 81 (2023) 102552 §4 に明記）。
Nd の評価はシンプロポーショネーション反応 (D) について +50.0 kJ/mol, +69.9 J/K/mol を
与え、これは Seifert/Rycerz の実測レンジ（+48〜55 kJ/mol, +63〜65 J/K/mol）と
整合するが、Pr の評価は +61.6 kJ/mol, +80.8 J/K/mol とレンジ外である。
**→ Nd の値の方が実測と整合しており、系列比較に使うなら Nd を基準にすべき。**

### 2-3. トレンド C：KLn₂Cl₇ は Ln が小さいほど出現しやすい

- La: KLa₂Cl₇ **存在しない**（代わりに K₃La₅Cl₁₈）
- Pr: KPr₂Cl₇ は KCl–PrCl₃ 系の評価に現れない（Cs 系には CsPr₂Cl₇ が存在）
- Nd: KNd₂Cl₇ 存在（742–783 K の狭い安定域）ΔH = +6.3 kJ/mol, ΔS = +39.9 J/K/mol
- Gd, Dy: KGd₂Cl₇, KDy₂Cl₇ とも存在
- Tb: KTb₂Cl₇ は一致融解する

すなわち **La 側では Ln-rich 側に K₃Ln₅Cl₁₈ 型（Ln の高配位数）が、
Er 側では KLn₂Cl₇ 型（低配位数）が現れる。** Seifert が Eu 系の抄録で
「Sm 化合物で現れた低配位数構造への傾向がさらに続く」と述べているのは
この点。**Deacon 条件で LnCl₃ が KCl とどの複塩を作るかは Ln により質的に変わる。**

### 2-4. Y vs Ho の決定実験に使えるデータ源

CLAUDE.md が挙げる Y（0.900 Å, 4f なし）と Ho（0.901 Å, 4f¹⁰）の対比には、
以下のペアがそのまま使える。**両者とも同じ Seifert グループ・同じ手法（DTA ＋
溶解熱量測定 ＋ 固体電解質 EMF）で測定されているため、系統誤差が相殺される。**

- Y: Seifert & Büchel, *Ternäre Chloride in den Systemen ACl/YCl₃ (A = Cs, Rb, K, Na)*,
  Z. Anorg. Allg. Chem. **624** (1998) 342–348.
- Ho: Roffe & Seifert, *Ternary chlorides in the systems ACl/HoCl₃ (A = Cs, Rb, K)*,
  J. Alloys Compd. **257** (1997) 128–133.
- 補助: Lu, Kang, He, Zhang, *Thermodynamic assessment of MCl–YCl₃ (M = Na, K, Rb, Cs) systems*,
  Calphad **47** (2014) 63–67（Y 系の CALPHAD パラメータ）。

### 2-5. 酸化還元活性 Ln（軸 A）について

- **Eu**: KCl 系で室温安定なのは K₂EuCl₅ のみ（K₃EuCl₆ は高温相、Cs/Rb 系にのみ
  T/H 型が確認されている）。EuCl₃ 自体が Deacon 条件で EuCl₂ に落ちるため、
  複塩データより先に Eu(III)/Eu(II) 平衡を見るべき。
- **Ce, Pr, Tb**: KCl–CeCl₃・KCl–PrCl₃ 系の複塩データは存在するが、酸化雰囲気では
  CeO₂/PrO₂ に行く（CLAUDE.md の前提どおり）。複塩データは「塩化物のまま
  留まった場合の下限シナリオ」としてのみ意味がある。

---

## 3. 出典ごとの詳細ノート

### 3-1. 【最重要・レビュー】Seifert の 2 部作レビュー（未入手）

- H. J. Seifert, *Ternary Chlorides of the Trivalent Early Lanthanides. Phase diagrams,
  crystal structures and thermodynamic properties*, **J. Therm. Anal. Cal. 67 (2002) 789–826**,
  DOI: 10.1023/A:1014341829611 （Ln = La–Gd）
- H. J. Seifert, *Ternary chlorides of the trivalent late lanthanides*,
  **J. Therm. Anal. Calorim. 83 (2006) 479–505**, DOI: 10.1007/s10973-005-7132-7 （Ln = Tb–Lu）

抄録から確認できた事実【抄録のみ】:
- 対象系: ACl/LnCl₃（A = Na, K, Rb, Cs）。K, Rb, Cs では A₃LnCl₆, A₂LnCl₅, ALn₂Cl₇。
  Ho 以降は Cs₃Ln₂Cl₉ も。Na では Na₃Ln₅Cl₁₈（La–Sm）, NaLnCl₄（Eu–Lu）。
- **手法**: ACl と LnCl₃ からの生成エンタルピーは **溶解熱量測定**。
  同反応の ΔG と ΔS は **温度可変の EMF 測定**（固体電解質ガルバニ電池）。
- 安定性の判定基準は「シンプロポーショネーション反応の自由エンタルピー ΔG°_syn < 0」。
- K₃LnCl₆（Ln = Ce, Pr, Nd）と Rb₃LaCl₆ は室温では存在せず、
  **ΔH = 48–55 kJ mol⁻¹ の格子エネルギー損失を ΔS = 63–65 J mol⁻¹ K⁻¹ の
  エントロピー利得が補償**して高温で生成する。
- Cs₃LnCl₆（Ln = La, Ce, Pr, Nd）は室温で安定で、固体内転移
  ΔH = 7.4–7.8 kJ mol⁻¹, ΔS = 10.9–11.2 J mol⁻¹ K⁻¹ を示す。

**PDF は Springer / akjournals いずれもペイウォールで取得できなかった（HTTP 403）。**
**本レビューを入手すれば、本レポートの「未取得」の大半が一挙に埋まる可能性が高い。**

### 3-2. 【全文確認】Sridar, Hao & Xiong (2023) — KCl–PrCl₃ / KCl–LiCl–PrCl₃

S. Sridar, L. Hao, W. Xiong, *Thermodynamic modeling of KCl-PrCl₃ and KCl-LiCl-PrCl₃ systems*,
**Calphad 81 (2023) 102552**, DOI: 10.1016/j.calphad.2023.102552.
全文 PDF: https://www.osti.gov/servlets/purl/2259214 （DOE PAGES, 無料）

Table 2（最適化パラメータ, J/mol）より **【全文確認】**:

```
G°(K2PrCl5) = 2 G°(KCl) + G°(PrCl3) − 50387
G°(K3PrCl6) = 3 G°(KCl) + G°(PrCl3) + 11199 − 80.774 T
```

Table 3（不変反応）は §1-1 の表に転記済み。

**重要な注意（§4 本文より）**: KCl-rich 側の不変反応
（`L → K₂PrCl₅ + K₃PrCl₆` および `K₃PrCl₆ → KCl + K₂PrCl₅`）が再現できなかったため、
**K₃PrCl₆ について文献の実験熱化学データ（refs 22–24）を評価に含めていない**。
同様に **K₂PrCl₅ の Seifert による実測生成エンタルピーも、PrCl₃-rich 側の不変反応
（`L → K₂PrCl₅ + PrCl₃`）が再現できないため除外している**。
→ **上の 2 式は「相図に合わせた値」であって「量熱測定値」ではない。**

引用されている一次熱化学文献:
- [20] H. J. Seifert, J. Sandrock, J. Uebach, *Zur Stabilität von Doppelchloriden in den
  Systemen ACl/PrCl₃ (A = Na–Cs)*, Z. Anorg. Allg. Chem. **555** (1987) 143–153,
  DOI: 10.1002/zaac.19875551215 ← K₂PrCl₅ の生成エンタルピー原典
- [21] M. Gaune-Escard, L. Rycerz, W. Szczepaniak, A. Bogacz, *Entropies of phase transitions
  in the M₃LnCl₆ compounds (M = K, Rb, Cs; Ln = La, Ce, Pr, Nd) and K₂LaCl₅*,
  J. Alloys Compd. **204** (1994) 189–192, DOI: 10.1016/0925-8388(94)90090-6
  ← **系列の転移エントロピーが直接載っている。優先入手候補。**
- [23] M. Gaune-Escard, L. Rycerz, *Heat Capacity of K₃LnCl₆ Compounds with Ln = La, Ce, Pr, Nd*,
  Z. Naturforsch. A **54** (1999) 229–235, DOI: 10.1515/zna-1999-3-412
- [24] T. Hattori, K. Igarashi, J. Mochinaga, *Enthalpies of Fusion of Intermediate Compounds,
  KMgCl₃, K₂MgCl₄, K₂BaCl₄, KCaCl₃, K₂SrCl₄, **K₂LaCl₅, K₃PrCl₆, K₃NdCl₆, KGd₃Cl₁₀, KDy₃Cl₁₀***,
  Bull. Chem. Soc. Jpn. **54** (1981) 1883–1884, DOI: 10.1246/bcsj.54.1883
  ← **融解エンタルピーの系列データ。BCSJ は J-STAGE で公開されている可能性が高い。**
  なお、この論文は Gd, Dy について **KGd₃Cl₁₀ / KDy₃Cl₁₀** という組成を挙げており、
  Seifert 系列の KLn₂Cl₇ と食い違う。要確認。

Z. Naturforsch. A 54 (1999) 229 の抄録（http://www.znaturforsch.com/aa/v54a/54a0229.pdf,
無料の 1 ページ抄録）は取得済みだが、**本文の数値表は含まれていない**。
zfn.mpdl.mpg.de（Vol.1 1946 – Vol.56 2001 が Open Access）に全文があるはず。

### 3-3. 【全文確認】Hao, Sridar & Xiong — KCl–LiCl–NaCl / KCl–LiCl–NdCl₃

L. Hao, S. Sridar, W. Xiong, *Thermodynamic description of molten salt systems:
KCl-LiCl-NaCl and KCl-LiCl-NdCl₃*, OSTI 2259217 / DOE PAGES.
全文 PDF: https://www.osti.gov/servlets/purl/2259217

Table 1（KCl-LiCl-NdCl₃ の最適化パラメータ, J/mol）より **【全文確認】**:

```
G°(K2NdCl5)      = 2 G°(KCl,solid) +   G°(NdCl3,solid) − 32900 − 18.57 T
G°(K3NdCl6)      = 3 G°(KCl,solid) +   G°(NdCl3,solid) + 17097 − 88.5  T
G°(KNd2Cl7)      =   G°(KCl,solid) + 2 G°(NdCl3,solid) +  6313 − 39.91 T
G°(LiNd2Cl7)     =   G°(LiCl,solid)+ 2 G°(NdCl3,solid) + 14251 − 19.831 T
G°(K3NdCl6,liq)  = 3 G°(KCl,liq)   +   G°(NdCl3,liq)   − 73900          （液相中の中性会合種）
```

液相モデル: `(K⁺, Li⁺, Nd³⁺)_P (Cl⁻, K₃NdCl₆)_Q`
— **中性会合種 K₃NdCl₆ を陰イオン副格子に入れて短距離秩序を表現している。**
これは Cu-K-Ln 融液を Temkin モデルより現実的に扱うときの直接の手本になる。

一次データの出典として引用されているのは:
- [34] Seifert, Fink, Uebach（KCl–NdCl₃ 状態図, DTA + XRD）
- [40] Gaune-Escard et al.（KCl–NdCl₃ 融液の混合エンタルピー, 1065 K）

### 3-4. 【全文確認】Gong et al. (2024/2025) — LiCl-KCl-LaCl₃ の Bayesian CALPHAD

R. Gong, S.-L. Shang, V. G. Goncharov, C. Cockrell, K. Trachenko, P. A. Bingham, X. Guo,
Z.-K. Liu, *Thermodynamic modeling of the LiCl-KCl-LaCl₃ system with Bayesian model
selection and uncertainty quantification*, arXiv:2406.15223v2.
PDF 取得・全文読了済み。

確認できた事実:
- **KCl–LaCl₃ 系の三元化合物は K₂LaCl₅ と K₃La₅Cl₁₈ のみ。K₃LaCl₆ も KLa₂Cl₇ も存在しない。**
- K₂LaCl₅ は Pnma（Meyer & Hüttl 1983）、K₃La₅Cl₁₈ は P6₃/m（Seifert）。
- 生成エンタルピーは Seifert et al. [26] の溶解熱量測定によるが、
  **本論文はその数値を引用文中に明示していない**（DFT + phonon で置き換えている）。
- 熱容量: K₂LaCl₅ の C_p は 500 K で 26.77 J/mol-K（Reuter & Seifert）、
  26.86 J/mol-K（Gaune-Escard & Rycerz）、26.99 J/mol-K（本論文 DFT-QHA）。
  K₃La₅Cl₁₈ は 500 K で 26.34 J/mol-K（Reuter & Seifert）、25.99（DFT-QHA）。
  ※ **単位は J/mol-atom-K の可能性が高い**（原子あたり）。要確認。
- 混合エンタルピー実測（Papatheodorou & Ostvold 1974）: 1173 K, x(LaCl₃)=0.496 で −15.319 kJ/mol。

一次データ出典:
- [26] H. J. Seifert, H. Fink, G. Thiel, *Thermodynamic properties of double chlorides in the
  systems ACl/LaCl₃ (A = Na, K, Rb, Cs)*, J. Less-Common Met. **110** (1985) 139–147,
  DOI: 10.1016/0022-5088(85)90315-7
- [27] G. Reuter, H. J. Seifert, *The heat capacities of ternary lanthanum chlorides AₙLaCl₃₊ₙ
  from 200 to 770 K and the ΔC_p values for their formation from nACl + LaCl₃*,
  Thermochim. Acta **237** (1994) 219–228, DOI: 10.1016/0040-6031(94)80178-9
  ← **ΔC_p が直接載っている唯一の文献。Neumann–Kopp 近似の妥当性検証に使える。**

### 3-5. 【抄録のみ】Seifert グループの Ln 別一次論文（ドイツ語抄録から確認）

いずれも **DTA で状態図 → 溶解熱量測定で ΔH（nACl + LnCl₃ 基準）→
固体電解質 EMF の温度依存性で ΔG, ΔS** という統一手法。抄録に定性的な結論は
書かれているが、**数値表は本文中で、いずれもペイウォール。**

| Ln | 文献 | 抄録から読める内容 |
|---|---|---|
| Pr | Seifert, Sandrock, Uebach, ZAAC **555** (1987) 143–153, DOI 10.1002/zaac.19875551215 | K₃PrCl₆, K₂PrCl₅ 等を確認。**室温で安定なのは A₂PrCl₅ のみ、他はすべて高温相**。La・Ce 化合物と同型 |
| Eu | Seifert & Sandrock, ZAAC **587** (1990) 110–118, DOI 10.1002/zaac.19905870113 | T/H-A₃EuCl₆（A = Cs, Rb）と T/H-NaEuCl₄ を新規発見。**KCl・NaCl 系で室温安定なのは A₂EuCl₅ と NaEuCl₄ のみ**。Sm 以降の低配位数構造への傾向が継続 |
| Gd | Seifert, Sandrock, Thiel, ZAAC **598** (1991) 307–318, DOI 10.1002/zaac.19915980128 | K, Rb, Cs 系は Sm³⁺・Eu³⁺ 系に対応。A₃GdCl₆, A₂GdCl₅, AGd₂Cl₇。K₂GdCl₅ は K₂PrCl₅ 型。溶解熱量＋EMF で隣接化合物からの ΔH°, ΔG° を決定 |
| Dy | Seifert & Krämer, ZAAC **620** (1994) 1543–1548, DOI 10.1002/zaac.19946200909 | A₃DyCl₆（エルパソライト族）と ADy₂Cl₇ が全アルカリで存在。A₂DyCl₅ は Cs, K のみ。**全化合物について (nACl + DyCl₃) からの生成エンタルピーを溶解熱量測定で決定**し、シンプロポーショネーションエンタルピーを算出。**K₃DyCl₆ のみが隣接相から格子エンタルピーを失って生成し、EMF = f(T) から 39 ℃ 以上でエントロピー安定化されることが判明** |
| Ho | Roffe & Seifert, J. Alloys Compd. **257** (1997) 128–133, DOI 10.1016/s0925-8388(96)03119-2 | 抄録取得できず |
| Y | Seifert & Büchel, ZAAC **624** (1998) 342–348, DOI 10.1002/(sici)1521-3749(199802)624:2<342::aid-zaac342>3.0.co;2-k | 抄録取得できず |
| Yb | Sebastian & Seifert, Thermochim. Acta **318** (1998) 29–37 | 抄録取得できず |
| Er | Dudek & Seifert, ZAAC **627** (2001) 2317–2322 | **K₃LnCl₆ は Ln = Ce–Lu および Y で存在（La は除く）**、その多形を整理 |
| Tb | Mitra, Uebach, Seifert, J. Solid State Chem. **115** (1995) 484–489, DOI 10.1006/jssc.1995.1163 | 抄録取得できず |
| Sm | **Thiel & Seifert, Thermochim. Acta 133 (1988) 275–282, DOI 10.1016/0040-6031(88)87169-7** | 抄録取得できず。**本プロジェクトの参照系** |
| Ce | Seifert, Sandrock, Thiel, J. Therm. Anal. **31** (1986) 1309–1318, DOI 10.1007/bf01914643 / Thiel, Sandrock, Seifert, Thermochim. Acta **92** (1985) 815–817 | 抄録取得できず |
| Nd | Seifert, Fink, Uebach, J. Therm. Anal. **33** (1988) 625–632, DOI 10.1007/bf02138565 | 抄録取得できず |

### 3-6. 【抄録のみ】Blachnik グループ（Seifert とは独立した一次データ源）

**Seifert とは別グループの独立測定であり、系統誤差の相互チェックに使える。**

- R. Blachnik, D. Selle, *Bildungsenthalpien von Alkalichlorid-Lanthanoidchlorid-Verbindungen*,
  Z. Anorg. Allg. Chem. **454** (1979) 82–89, DOI: 10.1002/zaac.19794540113
  > 抄録（独語）: 「A₃MCl₆ および AM₂Cl₇ 型（A = K, Cs; M = ランタノイド(III)）の一連の
  > 化合物の生成エンタルピーを、反応 `3 ACl + MCl₃ → A₃MCl₆` および
  > `ACl + 2 MCl₃ → AM₂Cl₇` に対して等周溶解熱量計で決定した。生成エンタルピーから
  > 円環過程を経て標準生成エンタルピーと格子エネルギーを算出し議論した。」
  → **反応定義が本レポートの (A), (C) と完全一致。系列の ΔH 表がここにあるはず。最優先。**

- R. Blachnik, D. Selle, *Zur Thermochemie von Alkalichlorid-Lanthanoid(III)-chloriden*,
  Z. Anorg. Allg. Chem. **454** (1979) 90–98, DOI: 10.1002/zaac.19794540114
  > 抄録: KCl–GdCl₃, KCl–DyCl₃（原文は「DyCl₂」だが誤植と思われる）, CsCl–PrCl₃/DyCl₃/ErCl₃/YbCl₃
  > の状態図を DTA で解明。A₃MCl₆, A₂MCl₅, AM₂Cl₇ の 3 種が存在。
  > 一致融解化合物の安定性の指標として「化合物の融点と外挿共晶温度の差」を提案。

- R. Blachnik, D. Selle, A. Schneider, *Die Bildungsenthalpie der Verbindung K₂LaCl₅*,
  Thermochim. Acta **3** (1971) 143–144, DOI: 10.1016/0040-6031(71)80008-4
  → **K₂LaCl₅ の生成エンタルピー単独論文。2 ページなので入手できれば即使える。**

- R. Blachnik, A. Jäger-Kasper, *Enthalpies of transformation of alkali halide–lanthanide(III)
  halide compounds of the type M₃LnX₆*, Thermochim. Acta **35** (1980) 259–262,
  DOI: 10.1016/0040-6031(80)87200-5
  → **M₃LnX₆ の転移エンタルピーの系列データ。トレンド A の定量化に直結。**

- F. Dienstbach, R. Blachnik, *Dampfdruckmessungen an Alkalichlorid-Gadoliniumchlorid-Schmelzen*,
  Z. Anorg. Allg. Chem. **442** (1978) 135–143, DOI: 10.1002/zaac.19784420117
  > 抄録: NaCl, KCl, CsCl と GdCl₃ の融液上の蒸気圧をトーションエフュージョン法で測定。
  > 凝縮生成物の分析から蒸気組成を決定し、**ACl, A₂Cl₂, GdCl₃, および複合気体種 AGdCl₄
  > の分圧を算出**（A = Na, K, Cs）。
  → **`volatility.py` に直結する。KLnCl₄(g) という複合気体種が Ln 側にも存在する
  > ことを示しており、Cu₃Cl₃(g) と同じ扱いが Ln についても必要かもしれない。**

### 3-7. 【抄録のみ】Rycerz & Gaune-Escard（相転移熱・熱容量）

- L. Rycerz, M. Gaune-Escard, *Enthalpies of Phase Transitions and Heat Capacity of TbCl₃
  and Compounds Formed in TbCl₃–MCl Systems (M = K, Rb, Cs)*,
  J. Therm. Anal. Calorim. **68** (2002) 973–981, DOI: 10.1023/a:1016102925181
  > 抄録: TbCl₃, KTb₂Cl₇, RbTb₂Cl₇, CsTb₂Cl₇, K₃TbCl₆, Rb₃TbCl₆, Cs₃TbCl₆ の
  > 固–固／固–液転移の molar enthalpy を DSC で決定。**M₃TbCl₆ と MTb₂Cl₇ は
  > いずれも一致融解し、加えて固–固転移を示す。** 固–固転移の Δ_trs H° は
  > M₃TbCl₆ で K: **6.1**、Rb: **7.6**、Cs: **7.0** kJ mol⁻¹。
  > MTb₂Cl₇ で Rb: **17.1**、Cs: **12.1 と 10.9** kJ mol⁻¹（K の値は抄録に記載なし）。
  > 熱容量は 300–1100 K で DSC 測定。
- M. Gaune-Escard, L. Rycerz, *Heat Capacity of K₃LnCl₆ Compounds with Ln = La, Ce, Pr, Nd*,
  Z. Naturforsch. A **54** (1999) 229–235（抄録ページのみ取得）
- L. Rycerz, M. Gaune-Escard, *Heat Capacity of the Rb₃LnCl₆ Compounds with Ln = La, Ce, Pr, Nd*,
  Z. Naturforsch. A **54** (1999) 397–403, DOI: 10.1515/zna-1999-6-709
- M. Gaune-Escard, L. Rycerz, *Unusual nature of lanthanide chloride–alkali metal chloride
  M₃LnCl₆ compounds in the solid state*, J. Alloys Compd. **408–412** (2006) 76–79,
  DOI: 10.1016/j.jallcom.2005.04.049

### 3-8. K–Cu 系

- **S. Niazi, A. Bonk, A. Hanke, M. to Baben, B. Reis, E. Olsen, et al.,
  *Thermal stability, hydrolysis and thermodynamic properties of molten KCl-CuCl*,
  Materialia 21 (2022) 101296, DOI: 10.1016/j.mtla.2021.101296.**
  Gold OA（Unpaywall で is_oa = true）だが、**ScienceDirect が全アクセス経路で 403 を返し、
  DLR elib のコピーは DLR 内部限定（401）で全文を取得できなかった。**
  抄録・検索スニペットから確認できた内容【抄録のみ】:
  - 40–80 mol% CuCl の組成域を冷却曲線と DSC で測定。
  - FactSage で CALPHAD 評価。液相は Redlich-Kister 多項式による subregular solution。
  - **共晶点 145.9 ℃、64.9 mol% CuCl。包晶温度 241.2 ℃。**
  - 中間化合物 **K₂CuCl₃** を化学量論相として扱い、**生成反応エンタルピーは
    ab initio 計算値で固定し、標準エントロピーを最適化**した。
  - 500 ℃ まで質量減少なし。加水分解試験で HCl 検出されず。
  - 生データ（DSC 等）は Zenodo に公開: https://zenodo.org/records/7287490
    （xlsx, 29.8 MB。**CALPHAD パラメータの実数値がここに入っている可能性がある**）

- **C. M. Fontana, E. Gorin, G. A. Kidder, C. S. Meredith,
  *Chlorination of Methane with Copper Chloride Melts. Ternary System CuCl–CuCl₂–KCl,
  and its Equilibrium Chlorine Pressures*, Ind. Eng. Chem. 44 (1952) 363–368,
  DOI: 10.1021/ie50506a044.**
  → **CLAUDE.md 作業優先順位 2 が求めている「CuCl₂-CuCl-KCl 融液の平衡塩素圧データ」
    そのもの。`fit_interactions()` の較正に直接使える。ACS ペイウォール（未入手）。**
  同シリーズ:
  - Fontana, Gorin, Kidder, Kinney, *…Oxygen Equilibrium Pressures and Oxide Solubility
    in the Melt*, Ind. Eng. Chem. 44 (1952) 369–373, DOI: 10.1021/ie50506a045
    → **酸素平衡圧と酸化物溶解度。Deacon 条件の p(O₂) 側にそのまま対応する。**
  - Gorin, Fontana, Kidder, Ind. Eng. Chem. 40 (1948) 2128–2134 / 2135–2138（速度・生成物分布）

- **H. Bloom, D. W. Williams, *A mass spectrometric study of the vapors above the molten salt
  systems LiCl–CuCl, KCl–CuCl, LiBr–CuBr, and NaI–CuI*, J. Chem. Phys. 75 (1981) 4636–4646,
  DOI: 10.1063/1.442579.**
  → **KCl–CuCl 融液上の気相種を質量分析で同定している。`data.py` の Cu₃Cl₃(g) の
    ΔHf/S° を検証する直接の実験データ源になりうる。要入手。**

- K. Żurowski, *Some aspects of the phase changes of the CuCl₂–KCl system*,
  J. Therm. Anal. 44 (1995) 197–204, DOI: 10.1007/bf02547148 【抄録のみ】
  > Cu/K モル比 0.2–2.0 の CuCl₂–KCl 混合物を空気および Ar 雰囲気下で TG 測定。
  > 試料は調製温度が異なるものを用いた。相変化と系の構成について結論。
- K. Żurowski, *Reproducibility of phase changes in the system CuCl₂–KCl*,
  J. Therm. Anal. 45 (1995) 437–445, DOI: 10.1007/bf02548776 【抄録のみ】
  > **533–543 K、573–593 K、603–623 K の相変化**の再現性を DTA ＋ X 線で調査。
  > 液体融液と固体状態について結論。
  → 653 K（380 ℃）はこれら全ての転移より上。**反応温度では CuCl₂–KCl 系は液相を含む。**

- S. Sutakshuto-Trivijitkasem, B. Holm, H. A. Øye, *The Phase Diagram CuCl₂–LiCl–KCl*,
  Acta Chem. Scand. **32a** (1978) 969–972, DOI: 10.3891/acta.chem.scand.32a-0969
  → Acta Chem. Scand. は全巻無料公開されているはずなので、**入手容易な唯一の
    Cu(II)–アルカリ塩化物状態図**。KCl–CuCl₂ 二元の境界情報が含まれる可能性が高い。

- **触媒側の直接文献（Deacon 反応そのもの）**:
  - A. Ya. Aglulin, *Effect of the Composition of Supported Copper-Containing Salt Catalysts
    on Their Activity in the Deacon Reaction: Dependence of the Rate of the Deacon Reaction
    on the Ratio between Copper and Potassium Chlorides in a Supported CuCl₂–KCl Salt Catalyst*,
    Kinet. Catal. **60** (2019) 290–296, DOI: 10.1134/s0023158419030017
  - A. Ya. Aglulin, *Kinetics and possible mechanism of HCl oxidation over supported
    copper-containing salt catalysts: II. Kinetics of HCl oxidation in the Deacon and methane
    oxychlorination reactions over the **CuCl₂–KCl–LaCl₃** catalyst*,
    Kinet. Catal. **55** (2014) 582–591, DOI: 10.1134/s0023158414050024
    → **Cu–K–La 系の速度論。本プロジェクトの前身系そのもの。**

- 無水 K₂CuCl₄ / KCuCl₃ について【二次資料経由（Wikipedia）】:
  - 無水 KCl/CuCl₂ 状態図では **KCuCl₃ が一致融解化合物として現れ、K₂CuCl₄ は現れない**。
  - 「無水 K₂CuCl₄ は複雑な経緯があり、存在しないかもしれない。加熱すると
    KCl·CuCl₂ と KCuCl₃ の混合物に変わる」
  - K₂CuCl₄·2H₂O は **93 ℃ 以上で KCl + KCuCl₃ + H₂O に分解**。
  - KCuCl₃（無水, 鉱物 belloite 型ではない方）は単斜晶 P2₁/c、
    a = 402.81, b = 1379.06, c = 873.35 pm, β = 97.137°, Z = 4、密度 2.86 g/cm³。
  - **これらの Wikipedia 記述の根拠となる一次文献（状態図論文）を特定できなかった。**

---

## 4. 要取得リスト（ペイウォール原典）

優先度順。★★★ = これが取れれば本レポートの空欄が最も埋まる。

| 優先 | 文献 | DOI | 何が取れるか |
|---|---|---|---|
| ★★★ | **Thiel G., Seifert H. J., *Properties of double chlorides in the systems ACl/SmCl₃ (A = Na–Cs)*, Thermochim. Acta 133 (1988) 275–282** | 10.1016/0040-6031(88)87169-7 | **参照系 Sm の K₃SmCl₆ / K₂SmCl₅ / KSm₂Cl₇ の ΔH, ΔG, ΔS と状態図** |
| ★★★ | **Blachnik R., Selle D., *Bildungsenthalpien von Alkalichlorid-Lanthanoidchlorid-Verbindungen*, ZAAC 454 (1979) 82–89** | 10.1002/zaac.19794540113 | **K₃MCl₆ と KM₂Cl₇ の系列 ΔH 表（Seifert と独立）** |
| ★★★ | **Seifert H. J., J. Therm. Anal. Cal. 67 (2002) 789–826（早期 Ln）** | 10.1023/A:1014341829611 | La–Gd の全複塩の ΔH, ΔG, ΔS と状態図の統合表 |
| ★★★ | **Seifert H. J., J. Therm. Anal. Calorim. 83 (2006) 479–505（後期 Ln）** | 10.1007/s10973-005-7132-7 | Tb–Lu, Y の同上 |
| ★★★ | **Fontana C. M. et al., Ind. Eng. Chem. 44 (1952) 363–368** | 10.1021/ie50506a044 | **CuCl–CuCl₂–KCl 融液の平衡塩素圧（`fit_interactions()` 較正用）** |
| ★★★ | **Fontana C. M. et al., Ind. Eng. Chem. 44 (1952) 369–373** | 10.1021/ie50506a045 | **同融液の平衡酸素圧・酸化物溶解度** |
| ★★ | Niazi S. et al., Materialia 21 (2022) 101296（Gold OA だが取得不可） | 10.1016/j.mtla.2021.101296 | K₂CuCl₃ の ΔH_f・S°、KCl–CuCl 状態図の全パラメータ |
| ★★ | Bloom H., Williams D. W., J. Chem. Phys. 75 (1981) 4636–4646 | 10.1063/1.442579 | **KCl–CuCl 融液上の気相種と分圧（Cu 揮発の一次データ）** |
| ★★ | Dienstbach F., Blachnik R., ZAAC 442 (1978) 135–143 | 10.1002/zaac.19784420117 | **KCl–GdCl₃ 融液上の p(ACl), p(A₂Cl₂), p(GdCl₃), p(AGdCl₄)** |
| ★★ | Gaune-Escard M. et al., *Entropies of phase transitions in the M₃LnCl₆ compounds (M = K, Rb, Cs; Ln = La, Ce, Pr, Nd) and K₂LaCl₅*, J. Alloys Compd. 204 (1994) 189–192 | 10.1016/0925-8388(94)90090-6 | **系列の転移エントロピー（トレンド A の定量化）** |
| ★★ | Blachnik R., Jäger-Kasper A., Thermochim. Acta 35 (1980) 259–262 | 10.1016/0040-6031(80)87200-5 | M₃LnX₆ の転移エンタルピー系列 |
| ★★ | Hattori T., Igarashi K., Mochinaga J., Bull. Chem. Soc. Jpn. 54 (1981) 1883–1884 | 10.1246/bcsj.54.1883 | K₂LaCl₅, K₃PrCl₆, K₃NdCl₆, KGd₃Cl₁₀, KDy₃Cl₁₀ の融解エンタルピー。**BCSJ は J-STAGE で公開の可能性大（未確認）** |
| ★★ | Seifert H. J., Krämer R., ZAAC 620 (1994) 1543–1548 | 10.1002/zaac.19946200909 | Dy 系の全 ΔH と 39 ℃ 転移の数値 |
| ★★ | Roffe M., Seifert H. J., J. Alloys Compd. 257 (1997) 128–133 | 10.1016/s0925-8388(96)03119-2 | **Ho 系（Y との対比実験の片割れ）** |
| ★★ | Seifert H. J., Büchel D., ZAAC 624 (1998) 342–348 | 10.1002/(sici)1521-3749(199802)624:2<342::aid-zaac342>3.0.co;2-k | **Y 系（Ho との対比実験の片割れ）** |
| ★ | Lu G., Kang Z., He M., Zhang Y., *Thermodynamic assessment of MCl–YCl₃ (M = Na, K, Rb, Cs) systems*, Calphad 47 (2014) 63–67 | 10.1016/j.calphad.2014.06.007 | K₃YCl₆, K₂YCl₅, KY₂Cl₇ の CALPHAD パラメータ |
| ★ | He M., Lu G., Kang Z., Zhang Y., *Thermodynamic assessment of the LiCl–KCl–CeCl₃ system*, Calphad 49 (2015) 1–7 | 10.1016/j.calphad.2015.01.006 | K₃CeCl₆, K₂CeCl₅, K₃Ce₅Cl₁₈ の CALPHAD パラメータ |
| ★ | Reuter G., Seifert H. J., Thermochim. Acta 237 (1994) 219–228 | 10.1016/0040-6031(94)80178-9 | **ΔC_p（Neumann–Kopp 近似の妥当性検証）** |
| ★ | Seifert H. J., Fink H., Thiel G., J. Less-Common Met. 110 (1985) 139–147 | 10.1016/0022-5088(85)90315-7 | K₂LaCl₅, K₃La₅Cl₁₈ の ΔH（溶解熱量測定） |
| ★ | Blachnik R., Selle D., Schneider A., Thermochim. Acta 3 (1971) 143–144 | 10.1016/0040-6031(71)80008-4 | K₂LaCl₅ の生成エンタルピー（2 ページ） |
| ★ | Seifert H. J., Sandrock J., Thiel G., ZAAC 598 (1991) 307–318 | 10.1002/zaac.19915980128 | Gd 系の全 ΔH, ΔG |
| ★ | Seifert H. J., Sandrock J., ZAAC 587 (1990) 110–118 | 10.1002/zaac.19905870113 | Eu 系の全 ΔH, ΔG |
| ★ | Sutakshuto-Trivijitkasem S., Holm B., Øye H. A., Acta Chem. Scand. 32a (1978) 969–972 | 10.3891/acta.chem.scand.32a-0969 | **CuCl₂–LiCl–KCl 状態図。Acta Chem. Scand. は無料公開のはずなので入手容易** |
| ★ | Aglulin A. Ya., Kinet. Catal. 55 (2014) 582–591 | 10.1134/s0023158414050024 | CuCl₂–KCl–LaCl₃ 触媒の速度論（前身系） |
| ★ | Kapała J. et al., *Prediction of the thermodynamic data of the pseudobinary alkali halide–lanthanide halide condensed systems*, J. Alloys Compd. 451 (2008) 679–681 | 10.1016/j.jallcom.2007.04.085 | **系列全体の推算 ΔH（記述子モデルの比較対象）** |
| ★ | Gaune-Escard M., Rycerz L., Ingier-Stocka E., Gadzuric S., *Compound formation in lanthanide–alkali metal halide systems*, Miner. Process. Extr. Metall. 123 (2013) 35–42 | 10.1179/0371955313z.00000000066 | 系列レビュー |

**入手できた無料全文（再取得不要）**
- Sridar, Hao & Xiong, Calphad 81 (2023) 102552 → https://www.osti.gov/servlets/purl/2259214
- Hao, Sridar & Xiong, KCl-LiCl-NaCl / KCl-LiCl-NdCl₃ → https://www.osti.gov/servlets/purl/2259217
- Gong et al., arXiv:2406.15223 → https://arxiv.org/pdf/2406.15223
- Gaune-Escard & Rycerz, Z. Naturforsch. A 54 (1999) 229 の抄録ページ →
  http://www.znaturforsch.com/aa/v54a/54a0229.pdf （**http のみ。https は証明書エラー**）
  Z. Naturforsch. の Vol.1–56 全文は https://zfn.mpdl.mpg.de/（Open Access）にあるはずだが、
  今回は個別記事の PDF URL パターンを特定できなかった。

---

## 5. 見つからなかったもの・ギャップ

### 5-1. 決定的に欠けているもの

1. **Sm 系の複塩データが皆無。** 参照系（Cu-K-Sm/γ-Al₂O₃, 25–30 wt% SmCl₃）に
   直接対応するデータが 1 件も数値として取れなかった。
   Thiel & Seifert (1988) が唯一の一次資料。
2. **K–Ln 複塩の ΔH を溶解熱量測定の一次値として取得できたものはゼロ。**
   取得できた 5 件（Pr 2 件、Nd 3 件）はすべて CALPHAD 最適化値であり、
   しかも Pr の評価は実測熱化学データを意図的に除外している。
   → **`data.py` に入れるなら `confidence` は低く設定し、`source` に
     「CALPHAD 最適化値、実測ではない」と明記すべき。**
3. **KCuCl₃ および K₂CuCl₄ の生成エンタルピー（KCl + CuCl₂ 基準）が
   実験値・評価値ともに見つからなかった。** K–Cu 側は KCl–CuCl（Cu(I)）の
   相図情報しか数値化できていない。Deacon 触媒で重要なのは Cu(II) 側なので、
   ここが最大の穴。
4. **KCuCl₃ の融点の数値。** 「一致融解する」という定性的記述しか得られなかった。
5. **無水 KCl–CuCl₂ 二元状態図の一次文献を特定できなかった。**
   Wikipedia 記述の出典を辿れていない。Sutakshuto-Trivijitkasem et al. (1978) が
   最有力候補。

### 5-2. データ間の食い違い（要解決）

- **KGd₃Cl₁₀ / KDy₃Cl₁₀ vs KGd₂Cl₇ / KDy₂Cl₇**:
  Hattori et al. (1981, BCSJ) は KGd₃Cl₁₀ と KDy₃Cl₁₀ の融解エンタルピーを報告しているが、
  Seifert 系列（ZAAC 598 (1991), ZAAC 620 (1994)）は KGd₂Cl₇, KDy₂Cl₇ を報告している。
  **どちらが正しいか、あるいは両方存在するか未確認。**
- **K₃NdCl₆ の室温安定性**: Seifert et al. は高温相（719 K 以上）とするが、
  Zhang et al. は室温まで安定で 2 つの多形転移を示すと報告。Ghosh et al. は
  急冷試料中に K₃NdCl₆ を検出できなかった。
  （OSTI 2259217 §2.3 に議論あり。同論文は Seifert を採用）
- **KNd₆Cl₁₉ という組成**が Hosoya et al. により報告されているが、
  XRD で確認されておらず、OSTI 2259217 の著者は否定的。
- **熱容量の単位**: arXiv:2406.15223 は K₂LaCl₅ の C_p を「26.99 J/mol-K」と
  書いているが、文脈（J/mol-atom-K を他所で使用）から **原子あたりの値**である
  可能性が高い。K₂LaCl₅ は 8 原子なので分子あたりなら約 216 J/mol/K。要確認。

### 5-3. 本調査でアクセスできなかった経路

- ScienceDirect（Elsevier）: すべて HTTP 403。Gold OA 論文でも同様。
- SpringerLink: 認証リダイレクトのみ。
- akjournals: **抄録ページは Python + ブラウザ UA で取得可能**（本文は不可）。
  → 他の JTAC 論文の抄録収集にはこの経路が使える。
- ResearchGate, Wiley 全文: 未試行または不可。
- Semantic Scholar API: HTTP 429（レート制限）。時間を空ければ使える可能性。
- Crossref API / OpenAlex API / Unpaywall API / OSTI API: **すべて正常動作**。
  特に **Crossref は ZAAC のドイツ語抄録を全文返す** ので、Seifert 系列の
  定性的結論の収集にきわめて有効だった。

---

## 6. `data.py` / `melt.py` への具体的な反映案

（本レポートの調査結果から示唆される事項。コードは変更していない。）

1. **K–Ln 複塩を `stability.py` の候補凝縮相に追加する場合**、380 ℃ では
   `K₂LnCl₅` のみが「実際に存在しうる固相」であり、`K₃LnCl₆` は
   Ln = La–Nd では安定域外（KCl + K₂LnCl₅ に分解）、Ln = Dy 以降では安定域内、
   という Ln 依存の場合分けが必要。境界は Sm–Gd 付近。
2. **`melt.py` の活量モデル**には、KCl–LnCl₃ 融液の短距離秩序を
   中性会合種 `K₃LnCl₆(l)` で表す方法が CALPHAD で実績がある
   （Nd 系で ΔG = −73.9 kJ/mol）。理想 Temkin からの第一の補正として
   Cu 側のクロロ銅酸錯体と並列に扱うとよい。
3. **`volatility.py`** には、Ln 側にも複合気体種 `KLnCl₄(g)` が存在する
   （Dienstbach & Blachnik 1978, Gd 系）ことを候補として入れるべき。
   KCl 添加が Cu 揮発を抑える一方で Ln 揮発を促進する可能性がある。
4. **`fit_interactions()` の較正**には Fontana et al. (1952) の
   CuCl–CuCl₂–KCl 平衡塩素圧が最適。同シリーズの酸素平衡圧論文と併せて
   使えば、`gas.py` の p(Cl₂)–p(O₂) と融液活量を一貫して結べる。
