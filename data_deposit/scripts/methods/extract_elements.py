#!/usr/bin/env python3
"""Extract element sets from xyz files under Original-Files/doi_files/.

For each DOI, walks leaf directories (those with no subdirs), picks ONE
representative *.xyz (preferring non-repacked), parses element symbols
from column 1, and produces:

  results/elements_per_paper.csv     doi, n_leaves, n_xyz_total, elements_seen
  results/element_histogram.csv      element, papers
  results/element_extraction_log.txt

Run:
  python3 extract_elements.py
"""

import os
import re
import sys
import csv
from collections import defaultdict, Counter
from pathlib import Path

DOI_ROOT = Path("/groups/bsavoie2/zli43/Original-Files/doi_files")
RESULTS  = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")
RESULTS.mkdir(exist_ok=True)

# Valid element symbols (H through Og, period 1-7). Used to filter junk
# tokens — some xyz files have headers or "X"/"TV" sentinel rows.
ELEMENTS = set("""
H He
Li Be B C N O F Ne
Na Mg Al Si P S Cl Ar
K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr
Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
""".split())

# Match an xyz coordinate line: <symbol> <x> <y> <z>
_NUM = r"-?\d+(\.\d+)?([eE][-+]?\d+)?"
_COORD_RE = re.compile(r"^\s*([A-Z][a-z]?)\s+({n})\s+({n})\s+({n})".format(n=_NUM))


def parse_elements(path, max_lines=2000):
    """Return set of element symbols parsed from xyz file.

    Reads up to max_lines (typical xyz is short, but a few are huge
    trajectory files; we just need the unique elements).
    """
    found = set()
    try:
        with open(path, errors="replace") as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    break
                m = _COORD_RE.match(line)
                if not m:
                    continue
                sym = m.group(1)
                if sym in ELEMENTS:
                    found.add(sym)
    except Exception:
        pass
    return found


def doi_from_path(p):
    parts = p.parts
    for i, part in enumerate(parts):
        if part.startswith("10.") and i + 1 < len(parts):
            return part + "/" + parts[i + 1]
        if part == "no_doi":
            return "no_doi/" + "/".join(parts[i+1:i+3])
    return None


def pick_one_xyz(leaf_dir):
    """Pick the first non-repacked xyz file in this leaf dir."""
    try:
        entries = sorted(leaf_dir.iterdir())
    except OSError:
        return None
    for e in entries:
        if e.name.endswith(".xyz") and not e.name.endswith("repacked.xyz") and e.is_file():
            return e
    return None


def main():
    # Step 1: find all leaf dirs (no subdirs, via -links 2 trick equivalent)
    # We use os.walk-once and detect leaves on the way.
    print("Walking doi_files to find leaf dirs ...", file=sys.stderr)
    leaves = []
    for dp, dnames, _ in os.walk(DOI_ROOT):
        if not dnames:
            leaves.append(Path(dp))
    print("  {} leaf dirs".format(len(leaves)), file=sys.stderr)

    # Step 2: for each leaf, pick one xyz, extract elements, attribute to DOI
    print("Parsing one xyz per leaf ...", file=sys.stderr)
    paper_elements = defaultdict(set)     # doi -> set(elements)
    paper_n_leaves = Counter()
    paper_n_xyz    = Counter()
    no_xyz_leaves  = 0
    parsed_leaves  = 0
    skipped_no_doi = 0

    for i, leaf in enumerate(leaves):
        if i % 20000 == 0:
            print("  {}/{} leaves ({} parsed)".format(
                i, len(leaves), parsed_leaves), file=sys.stderr)

        doi = doi_from_path(leaf)
        if doi is None:
            skipped_no_doi += 1
            continue

        # count total xyz files in this leaf for stats
        n_xyz_here = 0
        try:
            for e in leaf.iterdir():
                if e.name.endswith(".xyz") and not e.name.endswith("repacked.xyz"):
                    n_xyz_here += 1
        except OSError:
            pass
        paper_n_xyz[doi] += n_xyz_here
        if n_xyz_here == 0:
            no_xyz_leaves += 1
            continue

        sample = pick_one_xyz(leaf)
        if sample is None:
            no_xyz_leaves += 1
            continue

        els = parse_elements(sample)
        if els:
            paper_elements[doi] |= els
            parsed_leaves += 1
            paper_n_leaves[doi] += 1

    print("Done parsing.", file=sys.stderr)

    # Step 3: write per-paper element-set CSV
    out = RESULTS / "elements_per_paper.csv"
    with open(out, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["doi", "n_leaves_parsed", "n_xyz_total", "elements"])
        for doi in sorted(paper_elements):
            els = sorted(paper_elements[doi])
            w.writerow([doi, paper_n_leaves[doi], paper_n_xyz[doi], " ".join(els)])

    # Step 4: element frequency histogram (per paper)
    elem_papers = Counter()
    for doi, els in paper_elements.items():
        for el in els:
            elem_papers[el] += 1
    out2 = RESULTS / "element_histogram.csv"
    with open(out2, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["element", "papers"])
        for el, c in elem_papers.most_common():
            w.writerow([el, c])

    # Step 5: log
    with open(RESULTS / "element_extraction_log.txt", "w") as fp:
        fp.write("leaf dirs walked: {}\n".format(len(leaves)))
        fp.write("leaves parsed (had >=1 element): {}\n".format(parsed_leaves))
        fp.write("leaves with no xyz: {}\n".format(no_xyz_leaves))
        fp.write("leaves skipped (no DOI): {}\n".format(skipped_no_doi))
        fp.write("papers with at least one element: {}\n".format(len(paper_elements)))
        fp.write("unique elements seen: {}\n".format(len(elem_papers)))
        fp.write("\nTop 30 elements by paper count:\n")
        for el, c in elem_papers.most_common(30):
            fp.write("  {:<3s} {}\n".format(el, c))


if __name__ == "__main__":
    main()
