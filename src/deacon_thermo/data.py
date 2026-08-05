"""化学種のデータ。

数値の出典と信頼度は各エントリの `source` / `confidence` を見ること。
`DB.report()` で要検証のものが一覧できる。

一次資料:
  JANAF  = NIST-JANAF Thermochemical Tables, 4th ed.
  Barin  = I. Barin, Thermochemical Data of Pure Substances, 3rd ed.

LnOCl について
--------------
系列を通した実測値が揃わないため、既定では推定器を使う:

    dHf(LnOCl) = 1/3 dHf(Ln2O3) + 1/3 dHf(LnCl3) + DH_OXYCHLORIDE

Ln2O3 と LnCl3 は比較的よく分かっているので、この形にすると
**系列内の差**の精度が絶対値よりずっと良くなる（推定誤差が Ln 間で相関して
差分で相殺されるため）。単一の Ln の絶対的な安定性を論じるときは
この近似に頼らず一次資料を引くこと。
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
#: LnOCl 推定に使う反応エンタルピー [kJ/mol]
#:   1/3 Ln2O3 + 1/3 LnCl3 -> LnOCl
#: 系列を通して一定と仮定している。感度解析は sensitivity.py を参照。
#:
#: -58.0 は Sm の唯一の文献アンカー (Jacob 2016, Bull. Mater. Sci. 39, 603 の
#: Knacke 系 dG(1000 K) = -51.1 kJ/mol) を、下の dS = -8 J/mol/K と組で再現する値。
#: 実測系 4 元素 (La/Nd/Sm/Gd) でも同換算値は -35〜-60 とばらつくので、
#: 系列一定の仮定には ±15 kJ/mol 程度の元素依存誤差が残る
#: (docs/literature/lnocl_thermochemistry.md §G2)。
#: Koch & Cunningham (JACS 1952-54, 加水分解平衡の直接実測) の入手後に確定させる。
DH_OXYCHLORIDE = -58.0

#: LnOCl 推定に使う反応エントロピー [J/mol/K]（同反応）
#: Burns/Jacob 系の採用値 dS_av = -8 (±4)。旧仮定 +8 は符号が逆だった。
DS_OXYCHLORIDE = -8.0

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
        DB.add(
            Species(
                f"{ln}OCl(s)", "s",
                dHf_ox / 3 + dHf_cl3 / 3 + DH_OXYCHLORIDE,
                S_ox / 3 + S_cl3 / 3 + DS_OXYCHLORIDE,
                (72.0, 12.0, -5.0),
                {ln: 1, "O": 1, "Cl": 1},
                EST,
                f"推定: 1/3 {ln}2O3 + 1/3 {ln}Cl3 + {DH_OXYCHLORIDE} kJ/mol "
                f"(dS {DS_OXYCHLORIDE} J/mol/K)。絶対値は要検証、系列内の差にも"
                "±15 kJ/mol 程度の元素依存誤差がありうる (docs/literature 参照)",
            )
        )


_register_lanthanides()
