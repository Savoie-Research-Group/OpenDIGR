#!/usr/bin/env python3
"""Heatmap: xTB-IRC convergence rate, functional × element class.

Reads results/xtb_irc_convergence_functional_x_class.csv (which is
produced by xtb_irc_convergence.py) and writes:

  results/xtb_irc_convergence_heatmap.png

Cells: convergence rate (%) annotated with "rate%\n(n_conv / n_total)".
"""

import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

RESULTS = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")
SRC = RESULTS / "xtb_irc_convergence_functional_x_class.csv"

# Show only the named functionals with enough total N
MIN_TOTAL_N_FN = 1500
# Show only chemistry-classes (drop unclassified buckets)
KEEP_CLASSES = ["3d_TM", "4d_TM", "5d_TM",
                "lanthanides", "actinides", "main_group_only"]

# overall convergence rate (for the colormap centerpoint)
OVERALL = 0.1787


def main():
    rows = []
    with open(SRC) as fp:
        r = csv.DictReader(fp)
        for row in r:
            fn = row["functional"]
            cls = row["element_class"]
            n   = int(row["n_reactions"])
            cn  = int(row["n_converged"])
            rate = float(row["convergence_rate"])
            rows.append((fn, cls, n, cn, rate))

    # totals per functional (across kept classes only)
    fn_totals = defaultdict(int)
    fn_conv_totals = defaultdict(int)
    for fn, cls, n, cn, _ in rows:
        if cls not in KEEP_CLASSES:
            continue
        fn_totals[fn] += n
        fn_conv_totals[fn] += cn

    # candidate functionals: named and meeting MIN_TOTAL_N
    candidates = [(fn, fn_totals[fn]) for fn in fn_totals
                  if not fn.startswith("(no ") and fn_totals[fn] >= MIN_TOTAL_N_FN]
    # sort by N descending
    candidates.sort(key=lambda x: -x[1])
    # take top 15
    top_funcs = [c[0] for c in candidates[:15]]

    # also sort classes in a fixed chemistry-meaningful order
    classes = [c for c in KEEP_CLASSES]

    # build matrices
    rate_mat = np.full((len(top_funcs), len(classes)), np.nan)
    n_mat    = np.zeros_like(rate_mat, dtype=int)
    cn_mat   = np.zeros_like(rate_mat, dtype=int)
    lookup = {(fn, cls): (n, cn, rate) for fn, cls, n, cn, rate in rows}
    for i, fn in enumerate(top_funcs):
        for j, cls in enumerate(classes):
            if (fn, cls) in lookup:
                n, cn, rate = lookup[(fn, cls)]
                rate_mat[i, j] = rate * 100  # in %
                n_mat[i, j] = n
                cn_mat[i, j] = cn

    # Plot
    fig, ax = plt.subplots(figsize=(11, 8))
    # diverging cmap centered on overall rate (17.9%)
    cmap = plt.get_cmap("RdBu_r")
    # set color limits: cap at twice the overall rate so cells deep into red
    # stand out without being saturated; allow blue side same span.
    vmin = max(0, OVERALL*100 - 25)
    vmax = OVERALL*100 + 25
    im = ax.imshow(rate_mat, cmap=cmap, aspect="auto",
                   vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=20, ha="right")
    ax.set_yticks(range(len(top_funcs)))
    ax.set_yticklabels([f"{fn}  (n={fn_totals[fn]:,})" for fn in top_funcs])
    ax.set_xlabel("element class (chemistry-only, unclassified excluded)")
    ax.set_ylabel("DFT functional (paper's reported choice)")
    ax.set_title(
        f"xTB-IRC convergence rate (%), functional × element class\n"
        f"overall = {OVERALL*100:.1f}%; named functionals only with N≥{MIN_TOTAL_N_FN}; "
        f"cmap centered on overall rate"
    )

    # annotate cells
    for i in range(len(top_funcs)):
        for j in range(len(classes)):
            v = rate_mat[i, j]
            if np.isnan(v):
                continue
            n = n_mat[i, j]; cn = cn_mat[i, j]
            # text color: white when far from center
            text_color = "white" if abs(v - OVERALL*100) > 12 else "black"
            ax.text(j, i, f"{v:.1f}%\n({cn:,}/{n:,})",
                    ha="center", va="center", fontsize=7.5, color=text_color)

    cbar = fig.colorbar(im, ax=ax, label="convergence rate (%)")
    cbar.ax.axhline(OVERALL*100, color="black", linewidth=0.8)
    fig.tight_layout()
    out = RESULTS / "xtb_irc_convergence_heatmap.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
