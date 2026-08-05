# CuCl₂-CuCl-KCl(-LnCl₃) 融液上の平衡塩素圧 — 文献調査

調査日: 2026-08-05 / 調査手段: Web 検索・公開 PDF のみ（機関購読なし）

目的: `melt.py` の `fit_interactions()` に入れる `ClObservation`
（T, cu_total, diluents, p_Cl2, f_CuII）を作れる実測データを探す。

---

## 0. 結論（先に）

- **較正に直接使えそうな一次データセットは 4 件**見つかった（下表 A-D）。
  うち **KCl を含む系は A・B・C の 3 件**、**Ln（La）を含む系は B の 1 件のみ**。
- **最有力は B: Ruthven & Kenney, J. Inorg. Nucl. Chem. 30 (1968) 931-944。**
  理由は 3 つ:
  1. 測定量が p(Cl₂) そのもの（`kind="p"` の観測にそのまま入る）
  2. **CuCl₂-CuCl-KCl-LaCl₃ の四元系を含む唯一の実測**
  3. 抄録が「CuCl₂-CuCl と CuCl₂-CuCl-ZnCl₂ は本質的に理想、KCl を含む系は
     非理想」と明言しており、**本プロジェクトの「理想 Temkin が外れるのは
     クロロ銅酸錯体のせい」という仮説と測定事実の向きが一致している**
- ただし **A-D はいずれも本文 PDF を入手できていない**（ACS / Elsevier / Springer
  ペイウォール、WebFetch は 403）。**本報告に生の p(Cl₂) 数値は一つも載っていない。**
  数値は §6「要取得リスト」の手配が済むまで存在しないものとして扱うこと。
- 全文が取れて数値まで確認できたのは、**純 CuCl₂ の分解圧の「計算値」表（F）**と
  **溶融 KCl-CuCl₂ 触媒の相同定（G）**のみ。

---

## 1. 使えるデータセット一覧

| ID | 出典 | 系 | T | 組成範囲 | 測定量 | データ形式 | 信頼度 |
|----|------|-----|---|---------|--------|-----------|--------|
| **A** | Fontana, Gorin, Kidder, Meredith, *Ind. Eng. Chem.* **44** (1952) 363-368 | CuCl-CuCl₂-KCl（三元） | 未確認（メタン塩素化融液なので 400-500 ℃ 圏と推定、要確認） | 未確認 | 平衡塩素圧 + 相挙動（凝固点） | 未確認（表か図か不明） | 二次資料経由（抄録要約のみ） |
| **B** | Ruthven & Kenney, *J. Inorg. Nucl. Chem.* **30**(4) (1968) 931-944 | CuCl₂-CuCl / CuCl₂-CuCl-**KCl** / CuCl₂-CuCl-**KCl-LaCl₃** / CuCl₂-CuCl-ZnCl₂ | 未確認（同著者の速度論論文が 350-500 ℃ なので同程度と推定、要確認） | 未確認 | **平衡塩素圧**、そこから導いた混合エンタルピー・エントロピー | 未確認 | 二次資料経由（抄録要約のみ） |
| **C** | Shevelin, Molchanova, Yolshin, Batalov, *Electrochim. Acta* **48** (2003) 1385-1394 | CuCl-CuCl₂-**MeCl** (Me = Li, Na, **K**, Cs) | 800-1000 K | MeCl 0 → 50-90 mol%（Me 依存、K は約 60 mol% まで） | **Cu⁺ と Cu²⁺ の濃度**（+ 電子/イオン輸率、電導度、密度）、p(Cl₂) = 1 atm 固定 | 未確認 | 抄録のみ |
| **D** | Shevelin, Raskovalov, Molchanova, *Ionics* **23**(11) (2017) 3163-3168 | CuCl-CuCl₂（二元） | 833 K (560 ℃) | 二元のみ | **Cu⁺/Cu²⁺ 濃度比**を p(Cl₂) = 0.1-1 atm で | 未確認（本文中の図と推定） | 抄録のみ |
| E | *J. Solution Chem.* (2018), DOI 10.1007/s10953-018-0817-x「Physico-chemical Properties of the Molten CuCl-CuCl₂ System」 | CuCl-CuCl₂（二元） | 835, 866, 905, 943 K | 二元全域 | 密度・モル体積（最大泡圧法）、CuCl₂ 分解の平衡定数は**熱力学シミュレーション由来** | 抄録に数値あり（モル体積のみ） | 抄録のみ |
| F | Wang, Marin, Naterer, Gabriel「Thermodynamics and kinetics of the thermal decomposition of cupric chloride in its hydrolysis reaction」（Memorial Univ. リポジトリ PDF） | 純 CuCl₂(s) | 180-470 ℃ | 純物質 | ΔG と p(Cl₂) の関係（**NIST データからの計算値**）+ TGA 分解温度 | **表（全文確認済）** | 全文確認（ただし値は計算値） |
| G | Su, Mannini, Metiu, Gordon, McFarland, *Ind. Eng. Chem. Res.* **57** (2018) 7795-7801 | KCl-CuCl₂ 融液 (45:55 mol%) | 400, 450 ℃ | 1 点 | HCl 転化率、急冷後 XRD/SEM-EDS。**p(Cl₂) 測定なし** | 全文 PDF 確認済 | 全文確認 |

**除外したもの**: KCl-CuCl 二元の CALPHAD 評価（Niazi ら 2021, ScienceDirect PII
S2589152921002982）は Cu(I) のみで Cu(II) を含まず、塩素圧の情報を持たない。
LiCl-KCl-LaCl₃ の CALPHAD 評価群も Cu を含まないため対象外（ただし
W(LnCl₃-KCl) の目安としては後段で使える可能性あり）。

---

## 2. 抽出できた数値

### 2.1 純 CuCl₂ の分解: p(Cl₂) と自発分解温度（出典 F, Table 1）

CuCl₂(s) → CuCl(s) + ½Cl₂(g) の ΔG = 0 となる温度。

| p(Cl₂) [bar] | 10⁻⁶ | 10⁻⁴ | 10⁻² | 0.1 | 1 |
|---|---|---|---|---|---|
| 転移温度 [℃] | 180 | 260 | 360 | 425 | 470 |

- **これは測定値ではなく、NIST の Cu / Cl₂ / CuCl / CuCl₂ の熱力学量からの計算値。**
  出典 F の著者自身が「熱力学解析」と明記している。
- 同論文の TGA 実測分解開始温度は 410-420 ℃（Ar 掃気、20 ℃/min、
  試料 17-23 mg）。他者データとして固定床 390-450 ℃、DSC/TGA 450 ℃ を引用。
  著者は境界層効果で見かけ上 p(Cl₂) ≈ 1 bar に近い挙動になると説明している。
- **380 ℃ での対応する p(Cl₂)**: 上表の 360 ℃/10⁻² bar と 425 ℃/0.1 bar を
  log p 線形で内挿すると **約 0.02 bar**。ただしこれは筆者（本調査）が
  出典 F の表に対して行った内挿であって、出典 F にこの値は書かれていない。
  `melt.py` の `redox_K()` の桁の妥当性チェックにのみ使うこと。

### 2.2 溶融 KCl-CuCl₂ 触媒中の実在化学種（出典 G）

- 45 mol% KCl - 55 mol% CuCl₂、450 ℃、24 h 反応後に急冷した試料の XRD は
  **KCuCl₃、K₄Cu₄OCl₁₀（ponomarevite）、CuO** の 3 相。**CuCl は検出されず**。
- 反応前（新品）の急冷試料も同じ 3 相。反応後は K₄Cu₄OCl₁₀ と CuO が増加。
- 著者は KCuCl₃ と K₄Cu₄OCl₁₀ は冷却中に生成したものと解釈しているが、
  **急冷組織にクロロ銅酸塩（KCuCl₃）が出るという事実自体が、融液中で
  Cu-Cl-K の錯形成が起きていることの状況証拠**になる。
- 450 ℃、HCl:O₂ = 1:2、気泡塔高さ 8-11 cm で HCl 転化率 81%（平衡値 84%）。
  400 ℃ では 22%（液面接触のみ）→ 59%（10 cm）。24 h 活性低下なし。
- 反応 4（2HCl + ½O₂ = Cl₂ + H₂O）の ΔG は 400 ℃ で **-13.1 kJ/mol**（HSC 計算値）。

### 2.3 Cu 揮発の実測（出典: Catalysis in Industry 2016, DOI 10.1134/S2070050416040085）

- 350 ℃、5 時間での銅塩化物の蒸発損失:
  **Cu+K 塩化物系触媒で 0.45%**、**LaCl₃ 含有触媒で 0.72%**。
- 信頼度: **抄録／検索経由**。全文未確認。触媒担体・組成・気流条件が
  不明なので、この 2 数値の比較（La を入れると揮発が増える）を額面通りに
  受け取るのは危険。CLAUDE.md の「KCl/LnCl₃ は Cu 活量を下げる」という
  枠組みと逆向きに見えるので、**全文確認を優先すべき**。

### 2.4 Hisham & Benson, *J. Phys. Chem.* **99** (1995) 6194-6198（要検証）

OCR ミラーサイト経由で以下を得たが、**表番号・反応の対応付けが曖昧で、
そのまま `data.py` に入れてはいけない**。

- 塩素化ステップ（423 K）: ΔH° = -32.9 kcal/mol, ΔS° = -32.8 cal/(mol·K), ΔG° = -27.7 kcal/mol
- 酸化ステップ（673 K）: ΔH° = 26.1 kcal/mol, ΔS° = 4.7 cal/(mol·K), ΔG° = 22.9 kcal/mol
- 全反応 2HCl + ½O₂ = H₂O + Cl₂: ΔH° = -13.6 kcal, ΔS° = -15.4 cal/(mol·K)
- Table 3（実験観察）: CuO について「375 ℃ では Cl₂ の発生が遅い。
  420 ℃ では Cl₂ は速やかに発生するが触媒が気化した」

信頼度: **二次資料（OCR）経由・要検証**。ACS 原論文の Table 1 と突き合わせるまで
数値として採用しないこと。

### 2.5 工業条件（参考、US 4,119,705）

溶融銅塩化物 Deacon プロセスの特許請求範囲:
- KCl 20-40 wt%（好適 30 wt%）、残部が銅塩化物
- **希土類塩化物（好ましくは LaCl₃）を 5-20 wt%** 添加可
- 酸化反応器入口 820-870 °F（438-466 ℃）、4-6 atm
- 脱塩素反応器入口 930-970 °F（499-521 ℃）、1-1.5 atm
- 塩素収率 96-98%（HCl 基準）

平衡塩素圧の数値表は含まれていない。ただし **KCl 30 wt% + LaCl₃ 5-20 wt% という
組成が工業的に選ばれている**という事実は、参照系（25-30 wt% SmCl₃）の
妥当性の傍証になる。

---

## 3. 較正への適合性評価 — どれから手を付けるか

### 優先度 1: **B (Ruthven & Kenney 1968, JINC 30, 931)**

`melt.py` の設計と最も相性がよい。

- 測定量が p(Cl₂) なので `ClObservation(kind="p")` に直接入る。仕込み組成から
  f_CuII が決まる測定様式なら、`ClObservation` が想定している
  「T + 組成 + p(Cl₂)」の形にそのまま乗る。
- **KCl 系列と LaCl₃ 系列が同一著者・同一装置で取られている**ので、
  W(CuCl₂,KCl) と W(CuCl₂,LaCl₃) の**差**が系統誤差を打ち消して決まる。
  CLAUDE.md の「系列内の差は絶対値より信頼できる」原則にそのまま合う。
- ZnCl₂ 系が「理想」と報告されているのが効く。`fit_interactions()` の
  docstring が警告する縮退（W(CuCl₂,X) と W(CuCl,X) が差でしか決まらない）に対し、
  **「ZnCl₂ については両方 0」という参照点**が使えるので、KCl 系の非理想性が
  純粋に KCl 由来だと切り分けられる。
- 混合エンタルピー・エントロピーが論文中で既に導出されているので、
  正則溶液 W への換算がほぼ直読みでできる可能性がある（W ≈ ΔH_mix / y_i y_j）。
- **懸念**: 温度域。同著者の速度論論文が 350-500 ℃ なので 380 ℃ を含む可能性が
  高いが未確認。含まなければ外挿が必要で、W を温度非依存に置く現行仮定と
  組み合わせて誤差評価が要る。

### 優先度 2: **A (Fontana 1952, IEC 44, 363)**

- B と独立に KCl 三元の p(Cl₂) を持つ唯一のデータ。**B のクロスチェックとして必須。**
- 凝固点・溶解度データも同じ論文にあるので、液相線に進むとき（作業優先順位 4）に
  そのまま再利用できる。
- 相挙動の記述が本文にあるので、380 ℃ で液相が存在する組成域の判定にも使える。
- 姉妹論文（*Ind. Eng. Chem.* **44** (1952) 369-373、酸素平衡圧と融液中の酸化物溶解度）は
  **本プロジェクトの「オキシ塩化物か塩化物か」の問い**に直結する。優先度は同等。

### 優先度 3: **C (Shevelin 2003, Electrochim. Acta 48, 1385)**

- 唯一「**Cu²⁺ と Cu⁺ の濃度を独立に測っている**」データ。`ClObservation` の
  `kind="f"`（p(Cl₂) を条件に f_CuII を予測）に入る、稀な形式。
- p(Cl₂) = 1 atm 固定で **MeCl の種類（Li/Na/K/Cs）と濃度を振っている**ので、
  アルカリカチオンのサイズ依存が取れる。これは軸B（イオン半径）の考え方を
  カチオン側で検証する材料になる。
- **懸念が 2 つ**:
  1. 温度が 800-1000 K（527-727 ℃）で、380 ℃ から 150-350 ℃ 外挿になる。
  2. p(Cl₂) = 1 atm 固定なので、Cu(II)/Cu(I) 比が高い側に偏っている。
     実反応条件（p(Cl₂) ~ 0.01-0.1 atm 程度）とは Cu の酸化状態分布が違う。
- したがって **C は「W の符号と桁を独立に確認する」用途に留め、
  主フィットは B + A で行うのが妥当。**

### 優先度 4: **D (Ionics 2017)**

- 二元 CuCl-CuCl₂ のみなので W(CuCl₂,KCl) には効かない。
  だが **p(Cl₂) を 0.1-1 atm で振って Cu⁺/Cu²⁺ 比を測っている**ので、
  `redox_K(T)`（純液体基準の標準状態）そのものを 560 ℃ で検証できる。
  **F の計算値より価値が高いアンカー点。**

### 使わないほうがよいもの

- **E（J. Solution Chem. 2018）の「CuCl₂ 分解平衡定数」は熱力学シミュレーション由来**で
  実測ではない。密度・モル体積は実測だが本件には効かない。
- **F は計算値**。`redox_K()` の桁チェック以上には使えない。
- **G には p(Cl₂) の測定がない。** 相同定（KCuCl₃）は定性的証拠として使う。

### 縮退の扱い（実務メモ）

`fit_interactions()` の docstring どおり、塩素圧データは
a(CuCl₂)²/a(CuCl)² にしか感度がないので、KCl については
W(CuCl₂,KCl) と W(CuCl,KCl) の**差**しか決まらない。分離には

- KCl-CuCl 二元の混合熱（→ W(CuCl,KCl) を `fixed` で与える）
- または KCl-CuCl の CALPHAD 評価（Niazi ら 2021）の液相 Redlich-Kister 係数

が要る。後者は Cu(I) 側だけを扱っているので、まさに `fixed` に入れる相方として
適している。**これは B/A を入手する前でも着手できる。**

---

## 4. 出典詳細

### 一次データ（未入手・要手配）

1. **C. M. Fontana, E. Gorin, G. A. Kidder, C. S. Meredith**,
   "Chlorination of Methane with Copper Chloride Melts. Ternary System,
   CuCl-CuCl₂-KCl, and Its Equilibrium Chlorine Pressures",
   *Industrial & Engineering Chemistry* **44** (1952) 363-368.
   DOI: 10.1021/ie50506a044. 被引用 46（Semantic Scholar）。
   抄録は ACS により非公開。検索経由で得た内容要約:
   「CuCl-CuCl₂-KCl 三元融液の相挙動と融液上の平衡塩素圧を扱う。二元系
   CuCl-KCl および CuCl-CuCl₂ の凝固点・溶解度データ、および二元系
   CuCl-CuCl₂ 上の塩素圧は既報。」
   → **信頼度: 二次資料経由。**

2. **C. M. Fontana, E. Gorin, G. A. Kidder, R. Kinney**,
   "Chlorination of Methane with Copper Chloride Melts. Oxygen Equilibrium
   Pressures and Oxide Solubility in the Melt",
   *Industrial & Engineering Chemistry* **44** (1952) 369-373.
   DOI: 10.1021/ie50506a045. 被引用 12。
   → 融液中の酸化物溶解度と酸素平衡圧。**オキシ塩化物の議論に直結。**

3. **D. M. Ruthven, C. N. Kenney**,
   "Equilibrium Chlorine Pressures over Cupric Chloride Melts",
   *Journal of Inorganic and Nuclear Chemistry* **30**(4) (1968) 931-944.
   DOI: 10.1016/0022-1902(68)80312-4. 被引用 10（Semantic Scholar）。
   抄録は Elsevier により非公開。検索経由で得た内容要約:
   「亜酸化銅（cuprous）および第二銅（cupric）塩化物を含むハロゲン化物融液
   （液体系 CuCl₂-CuCl、CuCl₂-CuCl-KCl、**CuCl₂-CuCl-KCl-LaCl₃**、
   CuCl₂-CuCl-ZnCl₂）上の平衡塩素圧を測定した。データから混合エンタルピーと
   混合エントロピーを得た。CuCl₂-CuCl と CuCl₂-CuCl-ZnCl₂ は本質的に理想であり、
   塩化カリウムを含む混合物の非理想性の理由が議論されている。」
   → **信頼度: 二次資料経由（複数の検索結果で同一の要約が再現）。**

4. **D. M. Ruthven, C. N. Kenney**,
   "The Kinetics of the Oxidation of Hydrogen Chloride over Molten Salt Catalysts",
   *Chemical Engineering Science* **23**(9) (1968) 981-990.
   350-500 ℃、二元 CuCl-CuCl₂ / 三元 KCl-CuCl-CuCl₂ / **四元 KCl-LaCl₃-CuCl-CuCl₂**。
   K:Cu 等モル融液で見かけの活性化エネルギー **約 28 kcal/(mol HCl)**。
   **LaCl₃ 添加は酸素吸収ステップの触媒作用により活性を促進**すると結論。
   担持型 Deacon 触媒に結論が適用できると著者が明記。
   → 信頼度: 二次資料経由。

5. **D. M. Ruthven, C. N. Kenney**,
   "The Kinetics of Oxygen Absorption in Molten Salts Containing Cuprous Chloride",
   *Chemical Engineering Science* **22**(12) (1967) 1561-1570.
   → 3, 4 の前段。

6. **P. Yu. Shevelin, N. G. Molchanova, A. N. Yolshin, N. N. Batalov**,
   "Electron Transfer in an Electron-Ion Molten Mixture of CuCl-CuCl₂-MeCl
   (Me = Li, Na, K, Cs)", *Electrochimica Acta* **48** (2003) 1385-1394.
   DOI: 10.1016/S0013-4686(03)00005-7.
   抄録（検索経由）:「CuCl-CuCl₂-MeCl (Me = Li, Na, K, Cs) 融液の電子・イオン
   輸率、電気伝導度、密度、**Cu⁺ および Cu²⁺ 濃度**を、融液上の塩素分圧を
   大気圧に等しくして 800-1000 K で測定した。イオンおよび電子輸率の分係数と
   融液のモル体積を決定した。電子輸率は CuCl-CuCl₂ 融液で最大（tc = 0.9）。
   アルカリ金属塩化物を CuCl-CuCl₂ 融液に加えると、電子伝導成分は
   LiCl 90、NaCl 65、**KCl 60**、CsCl 50 mol% まで保たれる。」
   → 信頼度: 抄録のみ。

7. **P. Yu. Shevelin, A. A. Raskovalov, N. G. Molchanova**,
   "An Electron Transfer in CuCl-CuCl₂ Melt at Different Cl₂ Partial Pressures",
   *Ionics* **23**(11) (2017) 3163-3168. DOI: 10.1007/s11581-017-2120-z.
   抄録（検索経由）:「CuCl-CuCl₂ 融液の電子・イオン輸率、伝導度、
   **Cu⁺/Cu²⁺ 濃度比**を、温度 560 ℃、塩素分圧 0.1-1 atm で測定した。
   電子輸率と全伝導度は塩素分圧とともに増大…全伝導度は pCl₂ が 0.1 → 1 atm で
   約 6 → 約 8 S/cm、電子輸率は同区間で 67% → 95%。」
   → 信頼度: 抄録のみ。**Cu⁺/Cu²⁺ 比そのものの数値は抄録に無く、図の読み取りが必要。**

### 全文確認済み

8. **S. Su, D. Mannini, H. Metiu, M. J. Gordon, E. W. McFarland**,
   "Chlorine Production by HCl Oxidation in a Molten Chloride Salt Catalyst",
   *Ind. Eng. Chem. Res.* **57**(23) (2018) 7795-7801.
   DOI: 10.1021/acs.iecr.8b01141.
   OA PDF: https://cara.berkeley.edu/wp-content/uploads/2018/07/Industrial-Engineering-Chemistry-Research-2018-Su.pdf
   → **信頼度: 全文確認。** §2.2 の数値の出所。

9. **Z. Wang, G. Marin, G. F. Naterer, K. S. Gabriel**,
   "Thermodynamics and Kinetics of the Thermal Decomposition of Cupric Chloride
   in its Hydrolysis Reaction".
   PDF: Memorial University リポジトリ
   (memorial.scholaris.ca, bitstream 3c1943b4-53d0-4a8d-8d24-9ec03212259e)
   → **信頼度: 全文確認（PDF）。ただし掲載誌・巻号・年が PDF から特定できず。要確認。**
   §2.1 の Table 1 の出所。

### 周辺・傍証

10. **C. N. Kenney**, "Molten Salt Catalysis of Gas Reactions",
    *Catalysis Reviews* **11**(1) (1975) 197-224.
    → 3-5 の著者自身によるレビュー。**A・B の数値が図表で再録されている可能性が高い**。
      原論文が入手困難な場合の第一代替。

11. **J. Villadsen, H. Livbjerg**, "Supported Liquid-Phase Catalysts",
    *Catalysis Reviews* **17**(1) (1978) 203-272.
    → CuCl₂ 系 Deacon 触媒の調製・組成・構造、担持液相の分布モデル、
      液相内の拡散と反応。**A・B の再録の可能性がある第二の代替。**

12. **F. Wattimena, W. M. H. Sachtler**, "Catalyst Research for the Shell
    Chlorine Process", *Stud. Surf. Sci. Catal.* **7** (1981) 816-827.
    → **Cu-ジジム(didymium: Nd/Pr 混合)-K 塩化物 / SiO₂**。350-365 ℃。
      365 ℃・SV 120 L HCl/kg/h で平衡転化率 77% にほぼ等しい HCl 転化率。
      「活性物質の揮発」を克服すべき既知触媒の欠点として挙げている。
      **参照系 Cu-K-Sm/γ-Al₂O₃ の直接の先行系。Ln 系列比較の実データ源。**

13. **US 4,119,705**（Production of chlorine, 溶融塩 Deacon）
    → §2.5。

14. **Catalysis in Industry** (2016), DOI 10.1134/S2070050416040085,
    "Stability of catalysts for the oxidative chlorination of methane"
    → §2.3。**信頼度: 抄録／検索経由。**

15. **M. W. M. Hisham, S. W. Benson**, "Thermochemistry of the Deacon Process",
    *J. Phys. Chem.* **99**(16) (1995) 6194-6198. DOI: 10.1021/j100016a065.
    → §2.4。**OCR ミラー経由・要検証。**

16. **H. A. Andreasen, A. Mahan, N. J. Bjerrum**, "Densities of Molten Potassium
    Chloride-Copper(II) Chloride Obtained by the Automated Float Method",
    *J. Chem. Eng. Data* **26**(2) (1981) 195-197.
    → KCl-CuCl₂ 融液の密度。等モルで 2.225 g/cm³（出典 8 経由）。
      Bjerrum グループは溶融塩中のクロロ錯体の分光研究で著名なので、
      **同グループの KCl-CuCl₂ 分光論文を辿る価値がある**（§6 参照）。

17. **T. Zhang, C. Troll, B. Rieger, J. Kintrup, O. F. K. Schlüter, R. Weber**,
    "Composition Optimization of Silica-Supported Copper(II) Chloride Substance
    for Phosgene Production", *Appl. Catal. A* **365**(1) (2009) 20-27.
    → 出典 8 が **KCl-CuCl₂ 二元状態図**の出所として引用（45 mol% KCl で液相線 < 365 ℃）。

18. **G. N. Papatheodorou, O. J. Kleppa**, "Enthalpies of Mixing in the Liquid
    Mixtures of the Alkali Chlorides with MnCl₂, FeCl₂ and CoCl₂",
    *J. Inorg. Nucl. Chem.* (1971), PII 0022-1902(71)80419-0. 810 ℃（一部 690 ℃）。
    → **CuCl₂-KCl の混合熱そのものは見つからなかった**が、同族の 3d 二価塩化物 +
      アルカリ塩化物の W の桁を与える類推材料。

19. **KCl-CuCl₂ 系の中間化合物**: 無水 KCl/CuCl₂ 状態図では **KCuCl₃ が調和溶融**
    化合物として現れる。K₂CuCl₄ の無水相としての存在は歴史的に争われた
    （二水和物は 93 ℃ 以上で KCl + KCuCl₃ + H₂O に分解、無水化の成功は
    T. J. Nolan ら 1975 と記載）。
    → **信頼度: 三次資料（HandWiki）。要一次確認。** ただし
      **調和溶融する KCuCl₃ の存在は W(CuCl₂,KCl) < 0 の強い定性的根拠。**

20. **K. Niazi, A. Bonk, M. to Baben, B. Reis, E. Olsen, H. Nygård**（著者順要確認）,
    "Thermal stability, hydrolysis and thermodynamic properties of molten KCl-CuCl"
    (2021), ScienceDirect PII S2589152921002982。
    → KCl-CuCl **二元**の CALPHAD 評価（液相は Redlich-Kister 準正則、
      中間化合物 K₂CuCl₃ を化学量論相として扱う）。500 ℃ まで質量減少なし。
      **W(CuCl,KCl) を `fixed` で与えるための相方として最適。**

---

## 5. 要取得リスト（完全書誌）

優先順に。**1-3 が取れれば較正は動く。**

| # | 文献 | DOI / 識別子 | 取りたいもの |
|---|------|-------------|-------------|
| 1 | Ruthven & Kenney, *J. Inorg. Nucl. Chem.* 30(4) (1968) 931-944 | 10.1016/0022-1902(68)80312-4 | **p(Cl₂) の数値表（特に KCl 系と KCl-LaCl₃ 系）、測定温度、組成の刻み、導出された ΔH_mix / ΔS_mix** |
| 2 | Fontana, Gorin, Kidder, Meredith, *Ind. Eng. Chem.* 44 (1952) 363-368 | 10.1021/ie50506a044 | **CuCl-CuCl₂-KCl の p(Cl₂) 表、三元相図、測定温度域** |
| 3 | Fontana, Gorin, Kidder, Kinney, *Ind. Eng. Chem.* 44 (1952) 369-373 | 10.1021/ie50506a045 | **酸素平衡圧、融液中の酸化物溶解度**（オキシ塩化物問題に直結） |
| 4 | Shevelin ら, *Electrochim. Acta* 48 (2003) 1385-1394 | 10.1016/S0013-4686(03)00005-7 | **Cu⁺/Cu²⁺ 濃度の組成依存（Me = K の系列）**、測定温度ごとの値 |
| 5 | Shevelin ら, *Ionics* 23(11) (2017) 3163-3168 | 10.1007/s11581-017-2120-z | 560 ℃ の **Cu⁺/Cu²⁺ 比 vs p(Cl₂)** の図（読み取りが要る可能性大） |
| 6 | Ruthven & Kenney, *Chem. Eng. Sci.* 23(9) (1968) 981-990 | 10.1016/0009-2509(68)87084-8 | 速度式、四元系（KCl-LaCl₃）の組成、温度域 |
| 7 | Kenney, *Catal. Rev.* 11(1) (1975) 197-224 | 10.1080/01614947508079985 | **1・2 の数値の再録があるか**（原論文の代替） |
| 8 | Villadsen & Livbjerg, *Catal. Rev.* 17(1) (1978) 203-272 | 10.1080/03602457808080882 | 同上、および担持融液の分布モデル |
| 9 | Wattimena & Sachtler, *Stud. Surf. Sci. Catal.* 7 (1981) 816-827 | 10.1016/S0167-2991(09)60695-9（PII S0167299108646959 で ScienceDirect 上） | **Cu-ジジム-K 系の組成・寿命・揮発データ**（Ln 系列の一次データ） |
| 10 | Hisham & Benson, *J. Phys. Chem.* 99 (1995) 6194-6198 | 10.1021/j100016a065 | **Table 1 の原値**（§2.4 の検証） |
| 11 | Niazi ら (2021), KCl-CuCl | ScienceDirect PII S2589152921002982（掲載誌・巻号・DOI 未確定） | **液相 Redlich-Kister 係数**（W(CuCl,KCl) の固定用）、正式書誌 |
| 12 | *Catalysis in Industry* (2016) | 10.1134/S2070050416040085 | **揮発損失 0.45% / 0.72% の測定条件**（担体・組成・気流） |
| 13 | Wang, Marin, Naterer, Gabriel（CuCl₂ 分解） | 掲載誌不明 | **正式書誌**（PDF 本体は入手済） |
| 14 | Bjerrum グループの KCl-CuCl₂ 融液分光 | 未特定 | クロロ銅酸錯体の分光学的証拠と生成定数（Andreasen/Mahan/Bjerrum 1981 の引用文献から辿る） |
| 15 | KCl-CuCl₂ 無水二元状態図の一次資料 | Zhang ら *Appl. Catal. A* 365 (2009) 20-27 が引用 | KCuCl₃ の調和溶融、液相線 |

**入手手段の想定**: 1-6, 9-12 は機関購読または ILL。7, 8 は Taylor & Francis。
2, 3 は 1952 年の ACS 誌なので大学図書館の製本体でしか取れない可能性がある
（ACS Legacy Archives 契約の有無を確認）。

---

## 6. ギャップ（見つからなかったもの・存在しないと考えてよいもの）

1. **Sm を含む融液の平衡塩素圧データは存在しない。**
   CuCl₂-CuCl-KCl-SmCl₃ はもちろん、SmCl₃ を含む Cu 塩化物融液の熱力学測定
   （相図・混合熱・活量）を一件も見つけられなかった。
   → **参照系そのものの融液データは無い。**

2. **希土類を含む Cu 塩化物融液の平衡塩素圧は、Ruthven & Kenney (1968) の
   LaCl₃ 四元系ただ 1 件。** Nd, Gd, Ce, Y, Ho など他の Ln については皆無。
   → CLAUDE.md の「軸B（イオン半径）」を融液活量側から検証する実測は
     **La の 1 点しかない**。Ln 系列の予測は当面、記述子ベースの外挿に頼るしかない。

3. **380 ℃ という反応温度そのものでの融液 p(Cl₂) 実測は見つからなかった。**
   見つかった実測はいずれも 350 ℃ 以上、多くは 450-730 ℃。
   380 ℃ 近傍を直接カバーする可能性があるのは Ruthven & Kenney の 350-500 ℃ のみ。
   → **W の温度非依存仮定（`fit_interactions` の Notes）は、この外挿距離では
     妥当性を別途評価する必要がある。**

4. **CuCl₂-KCl 融液の混合エンタルピーの直接測定（カロリメトリ）が見つからない。**
   Papatheodorou & Kleppa は Mn/Fe/Co はやっているが Cu(II) はやっていない
   （Cu(II) 塩化物は測定温度で分解するので、通常のカロリメトリが難しいことが
   理由と推測される）。
   → **W(CuCl₂,KCl) と W(CuCl,KCl) の縮退を解く独立データが無い。**
     KCl-CuCl 側（Cu(I)）の CALPHAD 評価を `fixed` に使うのが現実的な回避策。

5. **溶融 CuCl₂-KCl 系のクロロ銅酸錯体の「生成定数」は見つからなかった。**
   水溶液系（Cu²⁺ + Cl⁻ の K₁, K₂ など）は多数あるが、融液中の値は無い。
   融液側の証拠は間接的なもの（調和溶融する KCuCl₃ の存在、急冷 XRD での
   KCuCl₃ 検出、KCl 60 mol% まで電子伝導が保たれること）に留まる。
   → **W を「錯形成の代理」として扱う現行の枠組み（`RegularSolution` の
     W < 0 が安定化に対応）は、この意味で妥当だが、W を独立の分光データで
     裏付けることは当面できない。**

6. **CuCl₂(s) 単体の分解圧の「実測」データが見つかっていない。**
   得られたのは計算値（出典 F）と TGA の分解開始温度のみ。
   Fontana 1952 が「二元 CuCl-CuCl₂ 上の塩素圧は既報」と述べているので、
   その引用文献に 1950 年以前の実測があるはず。
   → **要取得リスト #2 を入手したら、その参考文献を辿ること。**

7. **Deacon 触媒の速度論文献で Cu「活量」を定量している報告は見つからなかった。**
   K/La 添加の効果は、活性・Cu(II) 割合・分散度・凝集抑制として報告されており
   （例: K は CuCl → CuCl₂ の再生を促進し反応中の CuCl₂ 濃度を上げる、
   La は K-Cu の偏析と粒子凝集を抑制する）、
   **熱力学的活量に換算された値は無い。**
   揮発量の定量は §2.3 の 1 件のみ（しかも条件不明）。

---

## 7. 次アクション（本調査からの提案）

1. **要取得 #1（Ruthven & Kenney 1968 JINC）を最優先で手配する。**
   これ 1 件で KCl 系と LaCl₃ 系の両方が埋まり、かつ ZnCl₂ という
   「理想の参照点」が付いてくる。
2. **#11（KCl-CuCl の CALPHAD）は先行して着手できる。**
   W(CuCl,KCl) を固定できれば、#1 が来た時点で縮退なしに W(CuCl₂,KCl) が決まる。
3. #1 を入れる前に、**`ClObservation` を `kind="p"` で 2-3 点だけ手で作って
   `fit_interactions()` が回ることを確認しておく**とよい。
   `tests/test_model.py::test_ideal_model_reproduces_observed_lifetime` の
   xfail が W 導入で何桁動くかを先に見ておけば、
   #1 が到着したときに数値の妥当性を即座に判断できる。
4. **数値を `data.py` に入れる際は、CLAUDE.md の規律どおり
   「コミット → 図の再生成」の順を守ること。**
   本報告に載っている数値のうち `data.py` に入れてよいものは現時点で無い
   （§2.1 は計算値、§2.4 は要検証、§2.3 は条件不明）。
