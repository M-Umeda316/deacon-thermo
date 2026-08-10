"""Cu:K:Sm カチオン組成三角図(気相ポテンシャル拘束の相集合マップ)。

古典的な等温三元相図(液相線入り)ではない。CLAUDE.md のとおりこの系は
Cu-K-Sm-Cl-O-H 系なので、Cl と O を気相(380 C、HCl/O2 = 2:1 の平衡出口)の
ポテンシャルとして固定した「条件付き相図」を描く。色は Sm を宿す相の組で
塗り分ける(K-rich ほど Sm が塩化物複塩に係留される、が複塩仮説の絵)。

図は outputs/(git 管理外)に、データとの対応のためコミットハッシュ入りの
ファイル名で保存する(CLAUDE.md のデータ規律)。
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patheffects import withStroke
from matplotlib.tri import Triangulation

from deacon_thermo import cation_grid, gas_state

matplotlib.rcParams["font.family"] = ["Meiryo", "Yu Gothic", "sans-serif"]

#: Sm を宿す相の組 -> (凡例ラベル, 表示順)。順序は Sm-rich -> K-rich の空間順で、
#: 色は既定カテゴリパレットのスロット 1-5 を固定順で割り当てる(隣接ペア検証済み)。
CATEGORIES = [
    (frozenset({"SmOCl(s)"}), "SmOCl のみ(酸塩化物)"),
    (frozenset({"SmOCl(s)", "K2SmCl5(s)"}), "SmOCl + K₂SmCl₅ 共存"),
    (frozenset({"K2SmCl5(s)"}), "K₂SmCl₅ のみ(塩化物 100%)"),
    (frozenset({"K2SmCl5(s)", "K3SmCl6(s)"}), "K₂SmCl₅ + K₃SmCl₆"),
    (frozenset({"K3SmCl6(s)"}), "K₃SmCl₆ のみ(塩化物 100%)"),
]
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
GRAY = "#b8b8b8"  # Sm を含まない端(組成軸の縁)
INK = "#333333"


def category_index(assemblage) -> int:
    sm_hosts = frozenset(n for n in assemblage.phases if "Sm" in n)
    for i, (key, _) in enumerate(CATEGORIES):
        if sm_hosts == key:
            return i
    return len(CATEGORIES)  # Sm なし or 未分類


def to_xy(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """重心座標 (Cu, K, Sm) -> 平面。Cu 左下、K 右下、Sm 頂点。"""
    x = coords[:, 1] + 0.5 * coords[:, 2]
    y = (np.sqrt(3) / 2) * coords[:, 2]
    return x, y


def main() -> Path:
    gas = gas_state(653.15, hcl_o2_ratio=2.0)
    grid = cation_grid(gas, "Sm", n=80)
    cats = np.array([category_index(a) for a in grid.assemblages])
    x, y = to_xy(grid.coords)

    fig, ax = plt.subplots(figsize=(9.0, 8.2), dpi=150)
    cmap = ListedColormap(PALETTE + [GRAY])
    tri = Triangulation(x, y)
    # 三角形は頂点カテゴリの多数決で塗る(境界のにじみを抑える)
    tri_cat = np.array([np.bincount(cats[t], minlength=6).argmax() for t in tri.triangles])
    ax.tripcolor(tri, facecolors=tri_cat, cmap=cmap, vmin=-0.5, vmax=5.5)

    # 外枠と頂点ラベル
    ax.plot([0, 1, 0.5, 0], [0, 0, np.sqrt(3) / 2, 0], color=INK, lw=1.2)
    ax.text(-0.02, -0.03, "Cu", ha="right", va="top", fontsize=13, color=INK)
    ax.text(1.02, -0.03, "K", ha="left", va="top", fontsize=13, color=INK)
    ax.text(0.5, np.sqrt(3) / 2 + 0.025, "Sm", ha="center", va="bottom", fontsize=13, color=INK)

    # 主要領域の直接ラベル(色だけに頼らない)
    halo = [withStroke(linewidth=3, foreground="white")]
    for i, (_, label) in enumerate(CATEGORIES):
        mask = cats == i
        # 極小領域はラベルの重心が隣領域に食み出すので凡例に任せる
        if mask.sum() < 60:
            continue
        cx, cy = x[mask].mean(), y[mask].mean()
        short = label.split("(")[0].replace(" 共存", "")
        ax.text(cx, cy, short, ha="center", va="center", fontsize=9,
                color=INK, path_effects=halo)

    # 参照組成
    rx, ry = to_xy(np.array([[0.40, 0.35, 0.25]]))
    ax.plot(rx, ry, marker="*", ms=16, mfc="white", mec=INK, mew=1.5, ls="none")
    ax.annotate("参照組成 0.40:0.35:0.25", (rx[0], ry[0]), textcoords="offset points",
                xytext=(12, 10), fontsize=9, color=INK, path_effects=halo)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=c) for c in PALETTE] + [
        plt.Rectangle((0, 0), 1, 1, fc=GRAY)
    ]
    labels = [lab for _, lab in CATEGORIES] + ["Sm なし(組成軸の縁)"]
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.0, 1.0),
              fontsize=8.5, framealpha=0.9)

    commit = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True
    ).stdout.strip() or "unknown"
    ax.set_title(
        "Cu:K:Sm カチオン組成断面の安定相集合\n"
        f"380 ℃, HCl/O₂ = 2:1 平衡出口 (p(Cl₂) = {gas.p_Cl2:.3f} atm)。全域で CuCl₂(s) が共存",
        fontsize=11,
    )
    ax.text(0.99, -0.07,
            f"data.py @ {commit} / 固体候補のみ(Cu-K rich 側は実際には融体)",
            transform=ax.transAxes, ha="right", fontsize=7.5, color="#777777")
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()

    out = Path(__file__).resolve().parent.parent / "outputs"
    out.mkdir(exist_ok=True)
    path = out / f"fig_cation_triangle_sm_{commit}.png"
    fig.savefig(path, bbox_inches="tight")
    print(f"saved: {path}")
    return path


if __name__ == "__main__":
    main()
