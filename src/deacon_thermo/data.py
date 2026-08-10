"""化学種のデータ。

数値の出典と信頼度は各エントリの `source` / `confidence` を見ること。
`DB.report()` で要検証のものが一覧できる。

一次資料:
  JANAF  = NIST-JANAF Thermochemical Tables, 4th ed.
  Barin  = I. Barin, Thermochemical Data of Pure Substances, 3rd ed.

LnOCl について
--------------
LnOCl の G は反応量ベースの推定器で組む:

    dHf(LnOCl) = 1/3 dHf(Ln2O3) + 1/3 dHf(LnCl3) + DH_i
    S(LnOCl)   = 1/3 S(Ln2O3)   + 1/3 S(LnCl3)   + DS_i

(DH_i, DS_i) は元素別に文献値を持つ（LNOCL_PARAMS。Koch らの加水分解平衡実測、
Yang らの熱量測定）。文献値が無い元素は系列共通の既定値
(DH_OXYCHLORIDE, DS_OXYCHLORIDE) に落ちる。反応量（1/3 Ln2O3 + 1/3 LnCl3 →
LnOCl の ΔH, ΔS）で持つのは、文献ごとの補助データ系の違い（磁気エントロピーの
扱いなど）が反応量では相殺されるため。実測 DH_i は元素間で -35〜-60 kJ/mol と
大きくばらつくので、**系列一定の仮定は使わないこと**（旧設計の教訓）。
"""

from __future__ import annotations

from .species import Confidence, Database, Species

GOOD = Confidence.GOOD
FAIR = Confidence.FAIR
EST = Confidence.ESTIMATE

DB = Database()

# ---------------------------------------------------------------------------
# 気相
# ---------------------------------------------------------------------------
DB.add(Species("O2(g)", "g", 0.0, 205.15, (29.4, 3.5, -0.2), {"O": 2}, GOOD, "JANAF"))
DB.add(Species("Cl2(g)", "g", 0.0, 223.08, (36.9, 0.6, -2.8), {"Cl": 2}, GOOD, "JANAF"))
DB.add(
    Species(
        "HCl(g)", "g", -92.31, 186.90, (29.1, 1.0, 0.1), {"H": 1, "Cl": 1}, GOOD, "JANAF"
    )
)
DB.add(
    Species(
        "H2O(g)", "g", -241.83, 188.84, (30.0, 10.7, 0.3), {"H": 2, "O": 1}, GOOD, "JANAF"
    )
)

# 銅の気相種。CuCl 融体上では Cu3Cl3 が主要蒸気種だが、JANAF 値への差し替え後は
# 反応条件での支配種が CuCl2(g)（要検証のまま）に移っている点に注意。
DB.add(
    Species(
        "CuCl(g)", "g", 91.1, 237.0, (37.2, 0.3, -0.3), {"Cu": 1, "Cl": 1},
        GOOD, "JANAF 4th ed. (Chase 1998) table Cl-012 と一致確認 (91.086/237.207)",
    )
)
DB.add(
    Species(
        "Cu3Cl3(g)",
        "g",
        -258.571,
        429.526,
        (120.0, 1.0, -1.0),
        {"Cu": 3, "Cl": 3},
        GOOD,
        "JANAF 4th ed. (Chase 1998) table Cl-132。2026-08-05 に janaf.nist.gov で"
        "一次確認。Cp は暫定のまま (JANAF Cp298=124.57)。旧暫定値 -305/469 は"
        "純 CuCl 上の平衡圧が 100 atm になる誤りだった",
    )
)
DB.add(
    Species(
        "CuCl2(g)",
        "g",
        -43.0,
        278.0,
        (57.0, 0.5, -0.5),
        {"Cu": 1, "Cl": 2},
        EST,
        "要検証: JANAF/WebBook に気相データ無し。Cu3Cl3 差し替え後は反応条件での"
        "支配蒸気種なので揮発の結論を直接左右する。Wächter & Schäfer (1980) 要取得",
    )
)

# ---------------------------------------------------------------------------
# 銅の凝縮相
# ---------------------------------------------------------------------------
DB.add(
    Species(
        "CuCl(s)", "s", -137.2, 86.2, (38.3, 34.0, 0.0), {"Cu": 1, "Cl": 1},
        GOOD, "Barin", ((703.0, 10.2),),
    )
)
DB.add(
    Species(
        "CuCl2(s)", "s", -220.1, 108.1, (60.2, 34.0, 0.0), {"Cu": 1, "Cl": 2},
        GOOD, "Barin; 融点は分解が先行するため名目値", ((871.0, 20.0),),
    )
)
DB.add(
    Species("CuO(s)", "s", -157.3, 42.6, (48.6, 7.4, -7.6), {"Cu": 1, "O": 1}, GOOD, "JANAF")
)
DB.add(
    Species(
        "Cu2O(s)", "s", -168.6, 93.1, (62.3, 24.0, 0.0), {"Cu": 2, "O": 1}, GOOD, "JANAF"
    )
)
DB.add(
    Species(
        "Cu2OCl2(s)",
        "s",
        -389.0,
        145.0,
        (120.0, 20.0, 0.0),
        {"Cu": 2, "O": 1, "Cl": 2},
        EST,
        "melanothallite。操業点のすぐ近くに現れるので要検証",
    )
)

# ---------------------------------------------------------------------------
# アルカリ塩化物
# ---------------------------------------------------------------------------
DB.add(
    Species(
        "KCl(s)", "s", -436.5, 82.6, (41.4, 21.8, 3.2), {"K": 1, "Cl": 1},
        GOOD, "JANAF", ((1044.0, 26.3),),
    )
)

# ---------------------------------------------------------------------------
# ランタノイド
# ---------------------------------------------------------------------------
#: 反応 1/3 Ln2O3 + 1/3 LnCl3 -> LnOCl の系列共通の既定値。
#: **LNOCL_PARAMS に元素別値が無い元素にのみ使う**。感度解析（共通モードの
#: 不確かさ掃引）は sensitivity.py を参照。
#: -58.0 は Sm の Knacke 系アンカー (Jacob 2016, Bull. Mater. Sci. 39, 603 の
#: dG(1000 K) = -51.1 kJ/mol) を DS = -8 と組で再現する値で、Koch の直接実測
#: (653 K の dG と 0.1 kJ/mol 一致) により事後的に妥当性が確認された。
DH_OXYCHLORIDE = -58.0

#: 同反応の系列共通エントロピー既定値 [J/mol/K]。Burns/Jacob 系の dS_av = -8 (±4)。
DS_OXYCHLORIDE = -8.0

#: 元素別の LnOCl 反応量 (DH_i [kJ/mol], DS_i [J/mol/K] または None=既定値, 信頼度, 出典)。
#: 反応量で持つ理由は冒頭 docstring を参照。転記元の表はローカル管理
#: (docs/literature/extracted/、機密方針により git 管理外)。
LNOCL_PARAMS = {
    "La": (
        -60.4, +4.5, FAIR,
        "Koch, Broido & Cunningham, JACS 74 (1952) 2349 の平衡実測"
        "(JACS 75 (1953) 797 Table III に再掲) から当 DB の補助データで換算。"
        "Yang 2022 (Inorg. Chem. 61, 7590, Table 2, p.7593) の熱量測定 -73.5 とは"
        " 13 kJ/mol 不整合。反応温度域 (700-900 K) を直接測る Koch を優先",
    ),
    "Pr": (
        -63.5, +3.4, FAIR,
        "Koch & Cunningham, JACS 76 (1954) 1471-1474 要約 (p.1471) の"
        "加水分解平衡実測 (700-900 K) の 298 K 定数から換算。"
        "実測により dG_hyd(653 K) は既定推定の +20 から +6.6 kJ/mol に確定",
    ),
    "Nd": (
        -56.06, +2.47, FAIR,
        "Koch & Cunningham, JACS 76 (1954) 1471-1474 要約 (p.1471) から換算"
        "(La-Gd と手法統一)。独立な熱量測定 (Yang 2022 -60.49±0.42 +"
        " Gibson 2025 dS -7.1) と dG_hyd(653 K) で 2.3 kJ/mol 整合",
    ),
    "Sm": (
        -47.3, +7.5, FAIR,
        "Koch & Cunningham, JACS 75 (1953) 796-797, Table III, p.797 の"
        "加水分解平衡実測 (700-900 K) の 298 K 定数から換算。"
        "653 K の dG_hyd は実測と 0.1 kJ/mol で一致",
    ),
    "Gd": (
        -49.9, +3.4, FAIR,
        "同 Koch Table III, p.797。Yang 2022 の -48.57±0.70 と 1.3 kJ/mol で整合",
    ),
    "Ho": (
        -35.04, None, FAIR,
        "Yang 2022, Inorg. Chem. 61, 7590, Table 2, p.7593 (-35.04±0.70)。dS は既定値",
    ),
    "Er": (
        -37.83, None, FAIR,
        "Yang 2022, 同 Table 2 (-37.83±0.69)。dS は既定値",
    ),
    "Y": (
        -38.16, +9.4, FAIR,
        "dH: Yang 2022, 同 Table 2 (-38.16±0.77)。"
        "dS: Gibson 2025, J. Chem. Thermodyn. 211, 107549 (298 K 反応量)",
    ),
    # Ce, Eu, Dy は元素別実測なし -> 系列共通既定値 (EST)。
    # Ce は酸化雰囲気ではどのみち CeO2 に行く (軸A)。Eu/Dy は Burns 1998 の
    # 推定表 (0817、転記待ち) で置換予定。
}

#: Shannon イオン半径 (CN=6, 3+) [Å]
IONIC_RADIUS = {
    "La": 1.032, "Ce": 1.010, "Pr": 0.990, "Nd": 0.983, "Sm": 0.958,
    "Eu": 0.947, "Gd": 0.938, "Tb": 0.923, "Dy": 0.912, "Ho": 0.901,
    "Er": 0.890, "Tm": 0.880, "Yb": 0.868, "Lu": 0.861, "Y": 0.900,
}

#: 到達可能な酸化数（Deacon 条件下で機構に関わりうるもの）
#: 軸A: ここが空でない Ln は「第二の酸化還元中心」になりうる。
ACCESSIBLE_OXIDATION_STATES = {
    "Ce": (4,), "Pr": (4,), "Tb": (4,),
    "Eu": (2,), "Yb": (2,), "Sm": (2,),  # Sm(II) は酸化雰囲気では実質届かない
}

#: LnCl3 融点 [K]
LNCL3_MELTING_POINT = {
    "La": 1131, "Ce": 1090, "Pr": 1059, "Nd": 1031, "Sm": 955,
    "Eu": 896, "Gd": 882, "Tb": 855, "Dy": 927, "Ho": 991,
    "Er": 1049, "Y": 994,
}

# (dHf298 LnCl3, S298 LnCl3, dHf298 Ln2O3, S298 Ln2O3)  [kJ/mol, J/mol/K]
_LN_DATA = {
    "La": (-1071.1, 137.6, -1793.7, 127.3),
    "Ce": (-1053.5, 151.0, -1796.2, 148.1),
    "Pr": (-1058.6, 153.3, -1809.9, 152.7),
    "Nd": (-1041.8, 153.4, -1807.9, 158.6),
    "Sm": (-1025.9, 150.6, -1823.0, 150.6),
    "Eu": (-936.0, 144.1, -1651.4, 146.0),
    "Gd": (-1018.2, 151.4, -1819.7, 150.6),
    "Dy": (-993.1, 154.8, -1863.4, 149.8),
    "Ho": (-997.7, 158.2, -1880.7, 158.2),
    "Er": (-994.4, 154.8, -1897.9, 155.6),
    "Y": (-1000.0, 137.2, -1905.3, 99.1),
}

LANTHANIDES = tuple(_LN_DATA)


#: K-Ln 複塩の生成反応量。name -> (Ln, n_KCl, n_LnCl3, dHf [kJ/mol], dSf [J/mol/K],
#: transitions)。反応 n_KCl KCl(s) + n_LnCl3 LnCl3(s) -> 化合物(s)。
#: 出典: Seifert 2002, J. Therm. Anal. Cal. 67, 789-826, Table 8, p.817
#: (溶解熱量測定 + 固体電解質 EMF、dG は 573 K)。K0.6LnCl3.6 / K0.5LnCl3.5 の
#: LnCl3 規格化行は化学量論式に戻して 5 倍 / 2 倍してある。
#: K3LnCl6 (Ln = Sm, Eu, Gd) は低温形 (L) を基準に、L->H 転移を
#: T_tr = ddH/ddS で登録(表の L/H 両行から導出)。
#: K2CeCl5 / K2EuCl5 は原表の dG573 列が dH - T dS と不整合(誤植と判断し
#: dH/dS 列を採用。docs/literature/extracted/p2-1.csv の注記参照)。
#: K3LaCl6 と KLa2Cl7 は存在しない(表にも無い)ので登録しない。
K_LN_DOUBLE_SALTS = {
    "K2LaCl5(s)": ("La", 2, 1, -22.5, 15.1, ()),
    "K3La5Cl18(s)": ("La", 3, 5, -26.5, 47.0, ()),
    "K2CeCl5(s)": ("Ce", 2, 1, -29.1, 12.0, ()),
    "K3Ce5Cl18(s)": ("Ce", 3, 5, -24.5, 45.5, ()),
    "K3CeCl6(s)": ("Ce", 3, 1, 24.8, 79.9, ()),
    "K2PrCl5(s)": ("Pr", 2, 1, -31.4, 13.4, ()),
    "K3PrCl6(s)": ("Pr", 3, 1, 20.5, 82.2, ()),
    "K2NdCl5(s)": ("Nd", 2, 1, -35.0, 11.9, ()),
    "K3NdCl6(s)": ("Nd", 3, 1, 13.5, 79.8, ()),
    "KNd2Cl7(s)": ("Nd", 1, 2, 3.2, 34.0, ()),
    "K2SmCl5(s)": ("Sm", 2, 1, -41.8, 10.5, ()),
    "K3SmCl6(s)": ("Sm", 3, 1, -11.2, 60.7, ((606.1, 8.0),)),
    "KSm2Cl7(s)": ("Sm", 1, 2, -0.6, 42.0, ()),
    "K2EuCl5(s)": ("Eu", 2, 1, -44.3, 9.7, ()),
    "K3EuCl6(s)": ("Eu", 3, 1, -17.6, 58.7, ((630.4, 8.7),)),
    "KEu2Cl7(s)": ("Eu", 1, 2, -6.4, 41.4, ()),
    "K2GdCl5(s)": ("Gd", 2, 1, -49.7, 3.7, ()),
    "K3GdCl6(s)": ("Gd", 3, 1, -26.3, 51.3, ((635.7, 8.9),)),
    "KGd2Cl7(s)": ("Gd", 1, 2, -18.2, 29.8, ()),
}


def _register_double_salts() -> None:
    """複塩を二元塩化物からの合成で登録する。

    Cp は Neumann-Kopp(二元の和)なので、生成反応の dG(T) は厳密に
    dH - T dS になり、Seifert の EMF の表現(dG 線形)と整合する。
    二元側の融解転移は化合物には持ち込まない(いずれも固体で比較する)。
    """
    kcl = DB["KCl(s)"]
    for name, (ln, n_k, n_l, dh, ds, trans) in K_LN_DOUBLE_SALTS.items():
        lncl3 = DB[f"{ln}Cl3(s)"]
        DB.add(
            Species(
                name, "s",
                n_k * kcl.dHf298 + n_l * lncl3.dHf298 + dh,
                n_k * kcl.S298 + n_l * lncl3.S298 + ds,
                tuple(n_k * a + n_l * b for a, b in zip(kcl.cp, lncl3.cp, strict=True)),
                {"K": n_k, ln: n_l, "Cl": n_k + 3 * n_l},
                FAIR,
                f"{n_k} KCl + {n_l} {ln}Cl3 + ({dh} kJ/mol, {ds} J/mol/K)。"
                "Seifert 2002, JTAC 67, 789-826, Table 8, p.817 (EMF + 溶解熱量測定)",
                trans,
            )
        )
    # K-Cu(I) 複塩。380 C では融体側だが、相集合の候補としては登録しておく。
    cucl = DB["CuCl(s)"]
    DB.add(
        Species(
            "K2CuCl3(s)", "s",
            2 * kcl.dHf298 + cucl.dHf298 + (-8.09),
            2 * kcl.S298 + cucl.S298 + (-3.48),
            tuple(2 * a + b for a, b in zip(kcl.cp, cucl.cp, strict=True)),
            {"K": 2, "Cu": 1, "Cl": 3},
            EST,
            "2 KCl + CuCl + (-8.09 kJ/mol, -3.48 J/mol/K)。Niazi 2022, Materialia 21,"
            " 101296, Table 3, p.8 (dH は ab initio 固定、S は最適化) から換算。"
            "241 C 以上で分解(包晶)するので高温では融体の代理でしかない",
        )
    )


def _register_lanthanides() -> None:
    for ln, (dHf_cl3, S_cl3, dHf_ox, S_ox) in _LN_DATA.items():
        DB.add(
            Species(
                f"{ln}Cl3(s)", "s", dHf_cl3, S_cl3, (95.5, 13.4, 0.0),
                {ln: 1, "Cl": 3}, FAIR, "Barin / 系列内挿",
                ((float(LNCL3_MELTING_POINT.get(ln, 1000)), 48.0),),
            )
        )
        DB.add(
            Species(
                f"{ln}2O3(s)", "s", dHf_ox, S_ox, (114.0, 12.0, -14.0),
                {ln: 2, "O": 3}, FAIR, "Barin",
            )
        )
        dh, ds, conf, src = LNOCL_PARAMS.get(
            ln,
            (
                DH_OXYCHLORIDE, DS_OXYCHLORIDE, EST,
                "元素別実測なし。系列共通既定値 (±15 kJ/mol 程度の元素依存誤差あり)",
            ),
        )
        if ds is None:
            ds = DS_OXYCHLORIDE
        DB.add(
            Species(
                f"{ln}OCl(s)", "s",
                dHf_ox / 3 + dHf_cl3 / 3 + dh,
                S_ox / 3 + S_cl3 / 3 + ds,
                (72.0, 12.0, -5.0),
                {ln: 1, "O": 1, "Cl": 1},
                conf,
                f"1/3 {ln}2O3 + 1/3 {ln}Cl3 + ({dh} kJ/mol, {ds} J/mol/K)。{src}",
            )
        )


_register_lanthanides()
_register_double_salts()
