# 要取得文献 統合リスト(優先度順)

作成: 2026-08-05。`docs/literature/` の 4 報告書
(lnocl_thermochemistry / cl2_equilibrium_pressures / complex_salts_K_Ln_Cu /
cu_gas_species)の要取得リストを統合し、重複を除いて優先度を付け直したもの。
各文献の詳細な背景(抄録から分かっている内容、何に使うか)は元の報告書を参照。

優先度の基準:

- **S**: 結論(SmCl₃/SmOCl 判定、Cu 揮発の絶対値)を直接左右する、または較正を動かす
- **A**: 系列比較・複塩仮説・活量モデルの確定に必要
- **B**: 裏取り・補助データ・代替経路

---

## S: 最優先(これで結論が決まる)

| # | 文献 | DOI | 取りたいもの | 入手性 |
|---|------|-----|-------------|--------|
| S1 | Koch, C. W.; Cunningham, B. B. *J. Am. Chem. Soc.* **75**(4) (1953) 796–797 | 10.1021/ja01100a010 | **SmCl₃/GdCl₃ + H₂O = LnOCl + 2HCl の ΔH, ΔG, K(T)**。加水分解平衡の直接実測。DH_OXYCHLORIDE(現在 EST の −58)を確定させる決定打 | ACS Legacy Archives |
| S2 | Ruthven, D. M.; Kenney, C. N. *J. Inorg. Nucl. Chem.* **30**(4) (1968) 931–944 | 10.1016/0022-1902(68)80312-4 | **CuCl₂-CuCl(-KCl)(-KCl-LaCl₃)(-ZnCl₂) 融液の平衡塩素圧の数値表**、測定温度・組成刻み、導出済み ΔH_mix/ΔS_mix。`fit_interactions()` 較正の本命。ZnCl₂ 系が理想の参照点 | Elsevier |
| S3 | Fontana, C. M.; Gorin, E.; Kidder, G. A.; Meredith, C. S. *Ind. Eng. Chem.* **44** (1952) 363–368 | 10.1021/ie50506a044 | **CuCl-CuCl₂-KCl 三元の p(Cl₂) 表と相図・凝固点**。S2 の独立クロスチェック+液相線データ | ACS Legacy Archives |
| S4 | Fontana, C. M.; Gorin, E.; Kidder, G. A.; Kinney, R. *Ind. Eng. Chem.* **44** (1952) 369–373 | 10.1021/ie50506a045 | **同融液の酸素平衡圧・酸化物溶解度**。オキシ塩化物問題の融液側 | ACS Legacy Archives |
| S5 | Thiel, G.; Seifert, H. J. *Thermochim. Acta* **133** (1988) 275–282 | 10.1016/0040-6031(88)87169-7 | **K₃SmCl₆ / K₂SmCl₅ / KSm₂Cl₇ の ΔH, ΔG, ΔS と状態図**。参照系 Sm の複塩データ唯一の一次資料。K₃SmCl₆ の安定境界(380 ℃ を跨ぐのが Sm-Gd 付近)の裁定 | Elsevier |
| S6 | Wächter, H.; Schäfer, H. *Z. Anorg. Allg. Chem.* (1980)(巻・頁未確定) | — | **CuCl₂(g) 関連の平衡データ**。Cu₃Cl₃ の JANAF 差し替え後の支配蒸気種 CuCl₂(g)(EST)の唯一の手がかり。まず完全書誌の確定から | Wiley |
| S7 | Koch; Broido; Cunningham. *JACS* **74**(9) (1952) 2349–2351(La)/ Koch; Cunningham. *JACS* **76**(6) (1954) 1471–1474(Pr/Nd) | 10.1021/ja01129a049 / 10.1021/ja01635a003 | S1 の姉妹編。**La が系列比較の基準点、Pr/Nd は Knacke 系と JCT 2025 系の食い違い(Nd で 28 kJ/mol)の裁定** | ACS Legacy Archives |

## A: 系列・複塩・活量モデルの確定に必要

| # | 文献 | DOI | 取りたいもの | 入手性 |
|---|------|-----|-------------|--------|
| A1 | Seifert, H. J. *J. Therm. Anal. Cal.* **67** (2002) 789–826(La–Gd)/ **83** (2006) 479–505(Tb–Lu, Y) | 10.1023/A:1014341829611 / 10.1007/s10973-005-7132-7 | **K-Ln 複塩系列の ΔH, ΔG, ΔS と状態図の統合レビュー 2 部作**。これで複塩の「未取得」の大半が埋まる。個別 Ln 論文(下記 B 群)の上位互換 | Springer |
| A2 | Blachnik, R.; Selle, D. *Z. Anorg. Allg. Chem.* **454** (1979) 82–89 | 10.1002/zaac.19794540113 | **K₃MCl₆ / KM₂Cl₇ の系列 ΔH 表(等周溶解熱量計)**。Seifert と独立な測定系統でのクロスチェック | Wiley |
| A3 | Yang, S.; Anderko, A.; Riman, R. E.; Navrotsky, A. *Inorg. Chem.* **61**(19) (2022) 7590–7596 | 10.1021/acs.inorgchem.2c00763 | **全 REOCl の ΔfH₂₉₈ 数値表**(高温融体溶解熱量測定)。「系列内の差は信頼できる」前提(G2 問題)の裁定。HoOCl 実測が入っている見込み(Y/Ho 決定実験の裏打ち) | ACS |
| A4 | Gibson, A. et al. *J. Chem. Thermodyn.* **211** (2025) 107549 | 10.1016/j.jct.2025.107549 | **NdOCl / YOCl / TmOCl の S°₂₉₈ と Cp(T) 実測**、ΔG 規格化の定義。現行の S 推定(DS=−8)を実測に置換できる | Elsevier。**SSRN プレプリント 10.2139/ssrn.5208841 が先に取れる可能性** |
| A5 | Burns, J. B.; Peterson, J. R.; Haire, R. G. *J. Alloys Compd.* **265** (1998) 146–152 | 10.1016/S0925-8388(97)00435-0 | EuOCl/GdOCl/LuOCl の実測 ΔfH と、**全 Ln への推定表+使用した補助データ**(現行 DH 換算値の補助データ依存を解消) | Elsevier |
| A6 | Knacke, O.; Kubaschewski, O.; Hesselmann, K. *Thermochemical Properties of Inorganic Substances*, 2nd ed., Springer (1991) | 書籍 | **LaOCl/NdOCl/SmOCl/GdOCl の ΔHf₂₉₈, S°₂₉₈, Cp(T) 生データ**(Jacob 2016 の 1000 K 値の源流)。頁番号まで記録すること | 図書館 |
| A7 | Niazi, S. et al. *Materialia* **21** (2022) 101296 | 10.1016/j.mtla.2021.101296 | **KCl-CuCl 二元 CALPHAD(液相 Redlich-Kister 係数、K₂CuCl₃)**。W(CuCl,KCl) を `fixed` に入れて較正の縮退を解く。Gold OA だが直接取得は 403 だった | **Zenodo record 7287490 に生データ公開(即入手可)** |
| A8 | Bloom, H.; Williams, D. W. *J. Chem. Phys.* **75** (1981) 4636–4646 | 10.1063/1.442579 | **KCl-CuCl 融液上の気相種の質量分析**。Cu₃Cl₃(g)/CuCl(g) 分配の実測検証と、KCl 添加下での揮発種変化 | AIP |
| A9 | Seifert, H. J.; Büchel, D. *ZAAC* **624** (1998) 342–348(Y)/ Roffe, M.; Seifert, H. J. *J. Alloys Compd.* **257** (1997) 128–133(Ho) | 10.1002/(sici)1521-3749(199802)624:2<342::aid-zaac342>3.0.co;2-k / 10.1016/s0925-8388(96)03119-2 | **Y/Ho 対照ペアの複塩データ**(同グループ・同手法なので系統誤差が相殺)。`radius_controls()` の決定実験の熱力学的裏打ち | Wiley / Elsevier |
| A10 | Gaune-Escard, M. et al. *J. Alloys Compd.* **204** (1994) 189–192 | 10.1016/0925-8388(94)90090-6 | M₃LnCl₆(M=K,Rb,Cs; Ln=La–Nd)と K₂LaCl₅ の**転移エントロピー系列**(K₃LnCl₆ 安定境界トレンドの定量化) | Elsevier |
| A11 | Wattimena, F.; Sachtler, W. M. H. *Stud. Surf. Sci. Catal.* **7** (1981) 816–827 | 10.1016/S0167-2991(09)60695-9 | **Cu-ジジム(Nd/Pr)-K/SiO₂ の組成・寿命・揮発データ**。参照系の直接の先行系、Ln 系列比較の実測 | Elsevier |
| A12 | Guido, M.; Balducci, G.; Gigli, G.; Spoliti, M. *J. Chem. Phys.* **55** (1971) 4566 | 10.1063/1.1676789(要確認) | CuCl 蒸発の Knudsen 質量分析。**JANAF Cu₃Cl₃ 値の独立検証**(ΔvapH°(640 K)≈141 kJ/mol という検索要約値の真偽確認) | AIP |

## B: 裏取り・補助・代替経路

| # | 文献 | DOI / 識別子 | 取りたいもの |
|---|------|-------------|-------------|
| B1 | Kenney, C. N. *Catal. Rev.* **11** (1975) 197–224 | 10.1080/01614947508079985 | S2/S3 の数値が図表で再録されている可能性。**原論文が取れない場合の第一代替** |
| B2 | Villadsen, J.; Livbjerg, H. *Catal. Rev.* **17** (1978) 203–272 | 10.1080/03602457808080882 | 同上の第二代替+担持液相の分布モデル |
| B3 | Ruthven & Kenney *Chem. Eng. Sci.* **23**(9) (1968) 981–990 | 10.1016/0009-2509(68)87084-8 | KCl-LaCl₃ 四元系の速度論(La 添加は酸素吸収を促進)。S2 の姉妹編 |
| B4 | Shevelin et al. *Electrochim. Acta* **48** (2003) 1385–1394 | 10.1016/S0013-4686(03)00005-7 | Cu⁺/Cu²⁺ 濃度の実測(Me=Li/Na/K/Cs)。W の符号と桁の独立確認(800–1000 K からの外挿注意) |
| B5 | Shevelin et al. *Ionics* **23**(11) (2017) 3163–3168 | 10.1007/s11581-017-2120-z | 560 ℃ の Cu⁺/Cu²⁺ 比 vs p(Cl₂)。`redox_K(T)` のアンカー(図の読み取り要) |
| B6 | Hisham, M. W. M.; Benson, S. W. *J. Phys. Chem.* **99** (1995) 6194–6198 | 10.1021/j100016a065 | Deacon 熱化学の Table 1 原値(OCR 経由値の検証) |
| B7 | *Catalysis in Industry* (2016) | 10.1134/S2070050416040085 | **揮発損失 0.45%(Cu+K)/0.72%(LaCl₃ 含有)の測定条件**。「LnCl₃ は揮発を抑える」枠組みと逆向きに見えるため要確認 |
| B8 | Aglulin, A. Ya. *Kinet. Catal.* **55** (2014) 582–591(+ I 報 10.1134/S0023158414050012、2019 年の続報 10.1134/s0023158419030017) | 10.1134/s0023158414050024 | CuCl₂-KCl-LaCl₃ 触媒の速度論(前身系)。LaCl₃ で速度 1 桁向上の条件 |
| B9 | Blachnik, R.; Jäger-Kasper, A. *Thermochim. Acta* **35** (1980) 259–262 | 10.1016/0040-6031(80)87200-5 | M₃LnX₆ の転移エンタルピー系列 |
| B10 | Hattori, T.; Igarashi, K.; Mochinaga, J. *Bull. Chem. Soc. Jpn.* **54** (1981) 1883–1884 | 10.1246/bcsj.54.1883 | K₂LaCl₅, K₃PrCl₆, K₃NdCl₆, KGd₃Cl₁₀, KDy₃Cl₁₀ の融解エンタルピー。**J-STAGE で無料の可能性大**。KGd₃Cl₁₀ vs KGd₂Cl₇ の組成食い違いの確認も |
| B11 | Seifert & Krämer *ZAAC* **620** (1994) 1543–1548 | 10.1002/zaac.19946200909 | Dy 系複塩の全 ΔH と「39 ℃ 転移」の数値(トレンドの重希土側アンカー) |
| B12 | Seifert 系列の個別 Ln 論文(Pr: ZAAC 555 (1987) 143 / Nd: JTA 33 (1988) 625 / Ce: JTA 31 (1986) 1309 / Eu: ZAAC 587 (1990) 110 / Gd: ZAAC 598 (1991) 307 / Er: ZAAC 627 (2001) 2317 / Tb: JSSC 115 (1995) 484 / Yb: TCA 318 (1998) 29) | 各 DOI は complex_salts 報告書 §3-5 | A1 のレビュー 2 部作が取れれば**原則不要**。取れない場合の個別代替 |
| B13 | Reuter, G.; Seifert, H. J. *Thermochim. Acta* **237** (1994) 219–228 | 10.1016/0040-6031(94)80178-9 | AₙLaCl₃₊ₙ の ΔCp 実測(Neumann-Kopp 近似の妥当性検証) |
| B14 | Dienstbach, F.; Blachnik, R. *ZAAC* **442** (1978) 135–143 | 10.1002/zaac.19784420117 | KCl-GdCl₃ 融液上の p(KGdCl₄(g)) 等。**Ln 側の複合気体種**(K 添加が Ln 揮発を促進する可能性) |
| B15 | Lu et al. *Calphad* **47** (2014) 63–67(Y)/ He et al. *Calphad* **49** (2015) 1–7(Ce) | 10.1016/j.calphad.2014.06.007 / 10.1016/j.calphad.2015.01.006 | K-Y / K-Ce 複塩の CALPHAD パラメータ |
| B16 | Sutakshuto-Trivijitkasem et al. *Acta Chem. Scand.* **32a** (1978) 969–972 | 10.3891/acta.chem.scand.32a-0969 | CuCl₂-LiCl-KCl 状態図。**Acta Chem. Scand. は無料公開のはずで入手容易**。KCl-CuCl₂ 二元境界の一次候補 |
| B17 | Seifert, Fink, Thiel *J. Less-Common Met.* **110** (1985) 139–147 / Blachnik, Selle, Schneider *TCA* **3** (1971) 143–144 | 10.1016/0022-5088(85)90315-7 / 10.1016/0040-6031(71)80008-4 | K₂LaCl₅ の溶解熱量測定(後者は 2 頁のみ) |
| B18 | Guido, Gigli, Balducci *J. Chem. Phys.* **57** (1972) 3731 | 10.1063/1.1678845(要確認) | CuCl/Cu₂Cl₂ 解離エネルギー(JANAF CuCl(g) の相互検証) |
| B19 | Barin 3rd ed. (1995) 実データ頁(p.610 前後) | ISBN 978-3-527-28531-0 | CuCl₂ の気相データ収録有無の確定 |
| B20 | Gurvich / IVTANTHERMO の Cu-Cl 評価 | — | JANAF と独立な評価系統でのクロスチェック |
| B21 | Terlingen et al. *ACS Catal.* **12** (2022) 5698 の SI | OA: PMC9087184 | LnOCl→LnCl₃ 塩素化 ΔG の元データ。**無料、即入手可** |
| B22 | *Min. Metall. Explor.* (2021) carbochlorination レビュー | 10.1007/s42461-021-00490-z | Ln-O-Cl の ΔG 表と 400 ℃ 安定領域(Kellogg 図の独立チェック)。OA 版 hdl.handle.net/11336/151697 |
| B23 | Kapała, J. et al. *J. Alloys Compd.* **451** (2008) 679–681 | 10.1016/j.jallcom.2007.04.085 | アルカリ-Ln ハロゲン化物の系列推算 ΔH(記述子モデルの比較対象) |
| B24 | Gaune-Escard, M.; Rycerz, L. "Heat Capacity of K₃LnCl₆ Compounds with Ln = La, Ce, Pr, Nd" *Z. Naturforsch. A* **54** (1999) 229–235 | 10.1515/zna-1999-3-412 | K₃LnCl₆ の Cp(複塩を data.py に登録する際の Maier-Kelley 係数用)。OA だが bot 対策のためブラウザで取得のこと(De Gruyter or zfn.mpdl.mpg.de)。姉妹編の Rb₃LnCl₆ 版(54 (1999) 397、10.1515/zna-1999-6-709)は当面不要 |
| B25 | Wang, Marin, Naterer, Gabriel(CuCl₂ 分解) | 掲載誌未特定 | 正式書誌の確定のみ(PDF は入手済み) |
| B26 | Bjerrum グループの KCl-CuCl₂ 融液分光 | 未特定 | クロロ銅酸錯体の分光的証拠と生成定数。Andreasen/Mahan/Bjerrum *J. Chem. Eng. Data* **26** (1981) 195 の引用文献から辿る |
| B27 | 無水 KCl-CuCl₂ 二元状態図の一次資料 | Zhang et al. *Appl. Catal. A* **365** (2009) 20–27 が引用 | KCuCl₃ の一致融解・液相線(B16 が最有力候補) |

---

## 無料枠の取得結果(2026-08-05 実施)

取得物は `docs/literature/pdf/`(git 管理外)に保存。

| 対象 | 結果 | 備考 |
|---|---|---|
| Niazi 生データ(A7) | **取得済み** `niazi2022_kcl-cucl_rawdata.xlsx` | 3 組成(KCl 32/34/36 mol%)の融点・TGA・加水分解の生曲線 17 シート。**CALPHAD 係数は本文のみ**(依然 403)。生データから液相線を自前フィットする道あり |
| Terlingen SI(B21) | **取得済み** `terlingen2022_si/` + 本文 `terlingen2022_acscatal.pdf` | Europe PMC 経由。本文 13 頁も入手 |
| Acta Chem. Scand. 1978(B16) | **取得済み** `sutakshuto1978_cucl2-licl-kcl.pdf` | actachemscand.ki.ku.dk 公式。スキャン画像 PDF(テキスト層なし、読取は目視/OCR) |
| Gibson SSRN プレプリント(A4) | **取得済み**(ユーザーがブラウザで取得) `gibson2025_reocl_heat_capacity_ssrn.pdf` | 検証済み。Nd の実測 (ΔH −60.5, ΔS −7.1) が採用値 (−58, −8) を支持、Y は (−38.2, +9.4) で系列一定仮定の破れを実証。詳細は lnocl_thermochemistry.md の追記 |
| Z. Naturforsch. A(B24) | **スクリプト不可** | MPDL(Anubis)・De Gruyter(202 キュー)とも bot 対策。**OA なのでブラウザなら可** |
| Hattori BCSJ 1981(B10) | **無料経路消滅** | BCSJ の OUP 移管(2024)で J-STAGE アーカイブが撤去され、academic.oup.com でペイウォール化。**有料枠へ格下げ** |
| Min. Metall. Explor. OA 版(B22) | **リポジトリ停止中** | ri.conicet.gov.ar がタイムアウト。Unpaywall 上も OA 所在はここのみ。**後日リトライ** |

## 入手済み(再取得不要)

- Jacob, Dixit, Rajput, *Bull. Mater. Sci.* **39** (2016) 603(LnOCl の 1000 K 評価表)— PDF 取得・検証済み
- Sridar, Hao & Xiong, *Calphad* **81** (2023) 102552(K-Pr CALPHAD)— OSTI 2259214
- Hao, Sridar & Xiong(K-Nd CALPHAD)— OSTI 2259217
- Gong et al.(K-La, Bayesian CALPHAD)— arXiv:2406.15223
- Su et al., *IECR* **57** (2018) 7795(溶融 KCl-CuCl₂ 触媒)— OA PDF
- Wang et al.(CuCl₂ 分解、書誌未確定)— リポジトリ PDF

## 記録の規律(入手したら)

測定温度域・不確かさ・使用補助データを含め、`data.py` の `source` には
ページ番号まで書く(CLAUDE.md)。数値確定は個別コミットで行い、
該当する報告書の「未取得」欄にも追記して対応づけること。
