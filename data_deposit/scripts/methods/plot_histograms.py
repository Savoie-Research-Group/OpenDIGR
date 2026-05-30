#!/usr/bin/env python3
"""Top-N bar charts from the canonicalized per-paper CSVs.

Uses matplotlib (no seaborn). Writes PNGs to results/.
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")


def read_csv(path):
    rows = []
    with open(path) as fp:
        r = csv.reader(fp)
        next(r)  # header
        for row in r:
            if len(row) == 2:
                rows.append((row[0], int(row[1])))
    return rows


def bar(rows, title, xlabel, outpath, color="#4477AA", topn=20):
    rows = rows[:topn]
    if not rows:
        print("  (skip {}: no data)".format(outpath))
        return
    labels = [r[0] for r in rows]
    counts = [r[1] for r in rows]
    fig, ax = plt.subplots(figsize=(10, 7))
    bars = ax.barh(range(len(labels)), counts, color=color, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    # value labels
    xmax = max(counts) if counts else 1
    for i, (b, c) in enumerate(zip(bars, counts)):
        ax.text(c + xmax*0.005, i, "{:,}".format(c),
                va="center", ha="left", fontsize=9)
    ax.set_xlim(0, xmax * 1.12)
    ax.grid(axis="x", alpha=0.3)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print("wrote {}".format(outpath))


def heatmap_combos(rows, topn_func=15, topn_basis=15, outpath=None):
    """rows = list of (functional, basis_set, papers). Plot top-N x top-N heatmap."""
    from collections import Counter, defaultdict
    func_totals = Counter()
    basis_totals = Counter()
    combo = {}
    for fn, bs, c in rows:
        func_totals[fn] += c
        basis_totals[bs] += c
        combo[(fn, bs)] = c
    top_funcs = [f for f, _ in func_totals.most_common(topn_func)]
    top_basis = [b for b, _ in basis_totals.most_common(topn_basis)]
    M = [[combo.get((fn, bs), 0) for bs in top_basis] for fn in top_funcs]

    fig, ax = plt.subplots(figsize=(14, 9))
    im = ax.imshow(M, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(top_basis)))
    ax.set_xticklabels(top_basis, rotation=45, ha="right")
    ax.set_yticks(range(len(top_funcs)))
    ax.set_yticklabels(top_funcs)
    ax.set_xlabel("basis set")
    ax.set_ylabel("functional")
    ax.set_title("Top {} functionals × top {} basis sets (papers)".format(
        topn_func, topn_basis))
    # annotate cells with counts
    for i in range(len(top_funcs)):
        for j in range(len(top_basis)):
            v = M[i][j]
            if v == 0:
                continue
            colour = "white" if v > max(map(max, M))*0.55 else "black"
            ax.text(j, i, "{:,}".format(v),
                    ha="center", va="center", fontsize=8, color=colour)
    fig.colorbar(im, ax=ax, label="# papers")
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print("wrote {}".format(outpath))


def main():
    # Single-dimension histograms
    bar(read_csv(RESULTS / "canon_functional_per_paper.csv"),
        "Top 20 functionals (DFT and beyond)",
        "papers using this functional",
        RESULTS / "functional_top20.png", color="#4477AA")

    bar(read_csv(RESULTS / "canon_basis_set_per_paper.csv"),
        "Top 20 basis sets",
        "papers using this basis set",
        RESULTS / "basis_set_top20.png", color="#117733")

    bar(read_csv(RESULTS / "canon_software_per_paper.csv"),
        "Top 20 quantum chemistry software",
        "papers using this software",
        RESULTS / "software_top20.png", color="#CC6677")

    bar(read_csv(RESULTS / "canon_method_family_per_paper.csv"),
        "Top method families",
        "papers", RESULTS / "method_family_top20.png", color="#332288")

    bar(read_csv(RESULTS / "canon_calculation_types_per_paper.csv"),
        "Top calculation types",
        "papers", RESULTS / "calculation_types_top20.png", color="#882255")

    # Heatmap of functional × basis_set combos
    rows = []
    with open(RESULTS / "functional_x_basis_per_paper.csv") as fp:
        r = csv.reader(fp); next(r)
        for row in r:
            if len(row) == 3:
                rows.append((row[0], row[1], int(row[2])))
    heatmap_combos(rows, topn_func=15, topn_basis=15,
                   outpath=RESULTS / "functional_x_basis_heatmap.png")


if __name__ == "__main__":
    main()
