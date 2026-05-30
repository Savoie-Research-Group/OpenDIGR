#!/usr/bin/env python3
"""Re-render the xTB-IRC convergence heatmap with larger fonts for the SI."""
import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- CONFIGURE FOR YOUR SYSTEM (or pass as env vars) ---
# SRC : path to xtb_irc_convergence_functional_x_class.csv
# OUT : output path for the heatmap PNG
import os
SRC = Path(os.environ.get(
    "SRC",
    Path(__file__).resolve().parents[2] / "convergence" / "xtb_irc_convergence_functional_x_class.csv"
))
OUT = Path(os.environ.get(
    "OUT",
    Path.cwd() / "xtb_irc_convergence_heatmap.png"
))

MIN_TOTAL_N_FN = 1500
KEEP_CLASSES = ["3d_TM", "4d_TM", "5d_TM",
                "lanthanides", "actinides", "main_group_only"]
OVERALL = 0.1787

# Font sizes — bumped for SI readability
FS_CELL  = 12   # was 7.5
FS_TICK  = 13   # was ~10
FS_LABEL = 15   # was ~10
FS_TITLE = 15   # was ~10
FS_CBAR  = 14   # was ~10

def main():
    rows = []
    with open(SRC) as fp:
        r = csv.DictReader(fp)
        for row in r:
            rows.append((row["functional"], row["element_class"],
                         int(row["n_reactions"]), int(row["n_converged"]),
                         float(row["convergence_rate"])))

    fn_totals = defaultdict(int)
    for fn, cls, n, cn, _ in rows:
        if cls in KEEP_CLASSES:
            fn_totals[fn] += n

    candidates = [(fn, fn_totals[fn]) for fn in fn_totals
                  if not fn.startswith("(no ") and fn_totals[fn] >= MIN_TOTAL_N_FN]
    candidates.sort(key=lambda x: -x[1])
    top_funcs = [c[0] for c in candidates[:15]]

    classes = list(KEEP_CLASSES)
    rate_mat = np.full((len(top_funcs), len(classes)), np.nan)
    n_mat    = np.zeros_like(rate_mat, dtype=int)
    cn_mat   = np.zeros_like(rate_mat, dtype=int)
    lookup = {(fn, cls): (n, cn, rate) for fn, cls, n, cn, rate in rows}
    for i, fn in enumerate(top_funcs):
        for j, cls in enumerate(classes):
            if (fn, cls) in lookup:
                n, cn, rate = lookup[(fn, cls)]
                rate_mat[i, j] = rate * 100
                n_mat[i, j] = n
                cn_mat[i, j] = cn

    # Slightly bigger canvas so bigger labels don't crowd
    fig, ax = plt.subplots(figsize=(14, 10), dpi=150)
    cmap = plt.get_cmap("RdBu_r")
    vmin = max(0, OVERALL*100 - 25)
    vmax = OVERALL*100 + 25
    im = ax.imshow(rate_mat, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)

    # Cleaner class labels for x-axis (replace underscores with spaces)
    class_labels = [c.replace("_", " ") for c in classes]
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(class_labels, rotation=20, ha="right", fontsize=FS_TICK)
    ax.set_yticks(range(len(top_funcs)))
    ax.set_yticklabels([f"{fn}  (n={fn_totals[fn]:,})" for fn in top_funcs], fontsize=FS_TICK)
    ax.set_xlabel("Element class", fontsize=FS_LABEL)
    ax.set_ylabel("DFT functional (source paper)", fontsize=FS_LABEL)
    ax.set_title(
        f"xTB-IRC convergence rate (%), functional × element class\n"
        f"overall = {OVERALL*100:.1f}%; named functionals only with N≥{MIN_TOTAL_N_FN}; "
        f"colormap centered on overall rate",
        fontsize=FS_TITLE
    )

    for i in range(len(top_funcs)):
        for j in range(len(classes)):
            v = rate_mat[i, j]
            if np.isnan(v):
                continue
            n = n_mat[i, j]; cn = cn_mat[i, j]
            text_color = "white" if abs(v - OVERALL*100) > 12 else "black"
            ax.text(j, i, f"{v:.1f}%\n({cn:,}/{n:,})",
                    ha="center", va="center", fontsize=FS_CELL, color=text_color)

    cbar = fig.colorbar(im, ax=ax, label="convergence rate (%)")
    cbar.ax.axhline(OVERALL*100, color="black", linewidth=0.8)
    cbar.ax.tick_params(labelsize=FS_CBAR-1)
    cbar.set_label("convergence rate (%)", fontsize=FS_CBAR)
    fig.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")
    print(f"  size: {OUT.stat().st_size/1024:.0f} KB")

if __name__ == "__main__":
    main()
