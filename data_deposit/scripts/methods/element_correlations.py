#!/usr/bin/env python3
"""Join per-paper element sets (elements_per_paper.csv) with per-paper
functional / basis / software from comp_detail.json files.

Outputs (all in results/):
  paper_methods.csv               doi, functionals, basis_sets, softwares
  element_x_functional.csv        rows = element, cols = functional, vals = papers
  element_x_basis.csv             rows = element, cols = basis,       vals = papers
  element_class_breakdown.csv     class, papers, top functionals/basis
  element_x_functional.png        heatmap
  element_x_basis.png             heatmap
  element_class_breakdown.png     stacked bar by class
"""

import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Re-use canonicalization from sibling script.
sys.path.insert(0, str(Path(__file__).parent))
from extract_methods import (
    FUNCTIONAL_CANON, BASIS_CANON,
    normalize_value, canonicalize, canonicalize_software,
)

COMP_ROOT = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/doi_files")
RESULTS   = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")

# Element classifications
TM_3D = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
TM_4D = set("Y Zr Nb Mo Tc Ru Rh Pd Ag Cd".split())
TM_5D = set("Hf Ta W Re Os Ir Pt Au Hg".split())
LANTHANIDES = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
ACTINIDES   = set("Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr".split())
ALKALI       = set("Li Na K Rb Cs Fr".split())
ALKALINE_E   = set("Be Mg Ca Sr Ba Ra".split())
HALOGENS     = set("F Cl Br I At".split())
NOBLE_GASES  = set("He Ne Ar Kr Xe Rn".split())
ORGANIC_CORE = set("H C N O".split())   # we always have these
HEAVY_MAIN   = set("Tl Pb Bi Po".split())
METALLOIDS   = set("B Si Ge As Sb Te".split())

CLASSES = [
    ("3d transition metals",    TM_3D),
    ("4d transition metals",    TM_4D),
    ("5d transition metals",    TM_5D),
    ("Lanthanides",             LANTHANIDES),
    ("Actinides",               ACTINIDES),
    ("Halogens (heavy)",        HALOGENS - {"F"}),  # Cl/Br/I/At
    ("Heavy main-group (Tl-Po)",HEAVY_MAIN),
    ("Alkali metals",           ALKALI),
    ("Alkaline earth",          ALKALINE_E - {"Be"}),  # exclude tiny Be
]


def load_paper_elements():
    """Return {doi -> set(elements)}."""
    out = {}
    with open(RESULTS / "elements_per_paper.csv") as fp:
        r = csv.reader(fp); next(r)
        for row in r:
            doi, _, _, els = row
            out[doi] = set(els.split())
    return out


def doi_from_path(p):
    parts = p.parts
    for i, part in enumerate(parts):
        if part.startswith("10.") and i + 1 < len(parts):
            return part + "/" + parts[i + 1]
    return None


def extract_paper_methods():
    """Walk comp_details JSONs; per DOI return:
       (canonical functionals, basis sets, softwares) as sets.
    """
    funcs = defaultdict(set)
    bases = defaultdict(set)
    softs = defaultdict(set)

    json_files = list(COMP_ROOT.rglob("*-comp_detail.json"))
    print("  walking {} JSONs ...".format(len(json_files)), file=sys.stderr)
    for i, jp in enumerate(json_files):
        if i % 5000 == 0:
            print("    {}/{}".format(i, len(json_files)), file=sys.stderr)
        try:
            obj = json.load(open(jp))
        except Exception:
            continue
        entries = obj if isinstance(obj, list) else [obj]
        doi = doi_from_path(jp)
        if not doi:
            continue
        for e in entries:
            if not isinstance(e, dict):
                continue
            for v in normalize_value(e.get("functional")):
                funcs[doi].add(canonicalize(v, FUNCTIONAL_CANON))
            for v in normalize_value(e.get("basis_set")):
                bases[doi].add(canonicalize(v, BASIS_CANON))
            for v in normalize_value(e.get("software")):
                softs[doi].add(canonicalize_software(v))
    return funcs, bases, softs


def main():
    print("Loading per-paper elements ...", file=sys.stderr)
    paper_elements = load_paper_elements()

    print("Re-walking comp_details for per-paper methods ...", file=sys.stderr)
    funcs, bases, softs = extract_paper_methods()

    # ── save per-paper methods CSV
    all_dois = set(paper_elements) | set(funcs) | set(bases) | set(softs)
    with open(RESULTS / "paper_methods.csv", "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["doi", "functionals", "basis_sets", "softwares", "elements"])
        for doi in sorted(all_dois):
            w.writerow([
                doi,
                " | ".join(sorted(funcs.get(doi, set()))),
                " | ".join(sorted(bases.get(doi, set()))),
                " | ".join(sorted(softs.get(doi, set()))),
                " ".join(sorted(paper_elements.get(doi, set()))),
            ])

    # ── element frequency
    elem_papers = Counter()
    for doi, els in paper_elements.items():
        for el in els:
            elem_papers[el] += 1

    # ── element × functional matrix (only papers with BOTH info)
    both = set(paper_elements) & set(funcs)
    print("papers with both elements and functional info: {}".format(len(both)),
          file=sys.stderr)

    elem_func_count = defaultdict(lambda: Counter())  # element -> func -> papers
    elem_basis_count = defaultdict(lambda: Counter())
    for doi in both:
        for el in paper_elements[doi]:
            for fn in funcs[doi]:
                if fn:
                    elem_func_count[el][fn] += 1
            for bs in bases.get(doi, set()):
                if bs:
                    elem_basis_count[el][bs] += 1

    # ── select top elements + top functionals for heatmap
    TOP_ELEMS = 25
    TOP_METHODS = 18
    top_elems = [e for e, _ in elem_papers.most_common(TOP_ELEMS)]
    func_totals = Counter()
    for el in top_elems:
        for fn, c in elem_func_count[el].items():
            func_totals[fn] += c
    top_funcs = [f for f, _ in func_totals.most_common(TOP_METHODS) if f]

    basis_totals = Counter()
    for el in top_elems:
        for bs, c in elem_basis_count[el].items():
            basis_totals[bs] += c
    top_bases = [b for b, _ in basis_totals.most_common(TOP_METHODS) if b]

    # ── write CSV: element × functional
    def write_matrix(path, rows, cols, get):
        with open(path, "w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(["element"] + list(cols))
            for r in rows:
                w.writerow([r] + [get(r, c) for c in cols])

    write_matrix(RESULTS / "element_x_functional.csv", top_elems, top_funcs,
                 lambda e, f: elem_func_count[e].get(f, 0))
    write_matrix(RESULTS / "element_x_basis.csv", top_elems, top_bases,
                 lambda e, b: elem_basis_count[e].get(b, 0))

    # ── element-class breakdown
    class_papers     = Counter()
    class_func_top   = defaultdict(Counter)
    class_basis_top  = defaultdict(Counter)

    for doi, els in paper_elements.items():
        for cls_name, cls_set in CLASSES:
            if els & cls_set:                       # paper has any from class
                class_papers[cls_name] += 1
                for fn in funcs.get(doi, set()):
                    if fn:
                        class_func_top[cls_name][fn] += 1
                for bs in bases.get(doi, set()):
                    if bs:
                        class_basis_top[cls_name][bs] += 1

    with open(RESULTS / "element_class_breakdown.csv", "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["element_class", "papers",
                    "top1_functional", "top1_fn_papers",
                    "top2_functional", "top2_fn_papers",
                    "top3_functional", "top3_fn_papers",
                    "top1_basis", "top1_bs_papers",
                    "top2_basis", "top2_bs_papers",
                    "top3_basis", "top3_bs_papers"])
        for cls_name, _ in CLASSES:
            row = [cls_name, class_papers[cls_name]]
            for fn, c in class_func_top[cls_name].most_common(3):
                row += [fn, c]
            while len(row) < 8:
                row += ["", 0]
            for bs, c in class_basis_top[cls_name].most_common(3):
                row += [bs, c]
            while len(row) < 14:
                row += ["", 0]
            w.writerow(row)

    # ── plot element × functional heatmap
    def heatmap(rows, cols, M, title, xlabel, ylabel, outpath, cmap):
        fig, ax = plt.subplots(figsize=(max(12, len(cols)*0.7),
                                        max(8,  len(rows)*0.35)))
        im = ax.imshow(M, cmap=cmap, aspect="auto")
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=45, ha="right")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title)
        vmax = M.max() if M.size and M.max() else 1
        for i in range(len(rows)):
            for j in range(len(cols)):
                v = M[i][j]
                if v == 0:
                    continue
                colour = "white" if v > vmax*0.55 else "black"
                ax.text(j, i, "{:,}".format(int(v)),
                        ha="center", va="center", fontsize=7, color=colour)
        fig.colorbar(im, ax=ax, label="# papers")
        fig.tight_layout()
        fig.savefig(outpath, dpi=130)
        plt.close(fig)
        print("wrote {}".format(outpath))

    Mf = np.array([[elem_func_count[e].get(f, 0) for f in top_funcs]
                   for e in top_elems])
    heatmap(top_elems, top_funcs, Mf,
            "Functional usage by element (top 25 elements × top 18 functionals)",
            "functional", "element",
            RESULTS / "element_x_functional.png", "Blues")

    Mb = np.array([[elem_basis_count[e].get(b, 0) for b in top_bases]
                   for e in top_elems])
    heatmap(top_elems, top_bases, Mb,
            "Basis-set usage by element (top 25 elements × top 18 basis sets)",
            "basis set", "element",
            RESULTS / "element_x_basis.png", "Greens")

    # ── element-class stacked bar (papers per class)
    fig, ax = plt.subplots(figsize=(10, 6))
    classes = [c for c, _ in CLASSES]
    counts = [class_papers[c] for c in classes]
    bars = ax.barh(range(len(classes)), counts,
                   color="#CC6677", edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels(classes)
    ax.invert_yaxis()
    ax.set_xlabel("# papers containing any element from this class")
    ax.set_title("Papers by element class (transition metals dominate)")
    xmax = max(counts) if counts else 1
    for i, c in enumerate(counts):
        ax.text(c + xmax*0.005, i, "{:,}".format(c), va="center", fontsize=9)
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    ax.set_xlim(0, xmax*1.12)
    fig.tight_layout()
    fig.savefig(RESULTS / "element_class_breakdown.png", dpi=130)
    plt.close(fig)
    print("wrote {}".format(RESULTS / "element_class_breakdown.png"))

    # element frequency bar
    top_freq = elem_papers.most_common(40)
    fig, ax = plt.subplots(figsize=(12, 8))
    labels = [e for e, _ in top_freq]
    counts = [c for _, c in top_freq]
    ax.barh(range(len(labels)), counts, color="#117733",
            edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("# papers containing this element")
    ax.set_title("Top 40 elements across DIGR corpus")
    xmax = max(counts) if counts else 1
    for i, c in enumerate(counts):
        ax.text(c + xmax*0.005, i, "{:,}".format(c), va="center", fontsize=8)
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    ax.set_xlim(0, xmax*1.12)
    fig.tight_layout()
    fig.savefig(RESULTS / "element_frequency_top40.png", dpi=130)
    plt.close(fig)
    print("wrote {}".format(RESULTS / "element_frequency_top40.png"))


if __name__ == "__main__":
    main()
