#!/usr/bin/env python3
"""
build_tm_os_matrix.py

Aggregate transition_metal_oxidation_states.csv into a Metal x OS matrix
suitable for draw_tm_os_radial_dials_fullcircle.py.

Dedup rule per (parent_dir, stem-minus-last-two-_-tokens) group:
  1. Prefer the variant with charge == 0.
  2. Otherwise pick a single variant with |charge| == 1, preferring -1 over +1.
  3. If no |charge| <= 1 variant exists, skip the group.
  Within the chosen charge, prefer the lowest multiplicity; final tiebreak
  is alphabetical on the basename.

Sum reactant + product OS atom counts of the picked variant into per-metal bins.
OS values outside [NEG_MIN, POS_MAX] are clamped to the terminal bins.

Output: tm_os_matrix_SUM.csv with header
  Metal,-3,-2,-1,0,1,2,3,4,5,6,Total
"""

import csv
import os
import re
import sys
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
IN_CSV = os.path.join(ROOT, "transition_metal_oxidation_states.csv")
OUT_CSV = os.path.join(HERE, "tm_os_matrix_SUM.csv")

NEG_MIN = -3
POS_MAX = 6
OS_COLS = list(range(NEG_MIN, POS_MAX + 1))

TM_ORDER = [
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
]
TM_SET = set(TM_ORDER)

ATOM_RE = re.compile(r"([A-Z][a-z]?)(\d+):(-?\d+)")


def parse_stem(basename: str):
    """Return (group_stem, charge, multiplicity) or None if malformed."""
    if not basename.endswith(".tar.zst"):
        return None
    stem = basename[:-len(".tar.zst")]
    parts = stem.rsplit("_", 2)
    if len(parts) != 3:
        return None
    try:
        charge = int(parts[1])
        mult = int(parts[2])
    except ValueError:
        return None
    return parts[0], charge, mult


def clamp(k: int) -> int:
    if k <= NEG_MIN:
        return NEG_MIN
    if k >= POS_MAX:
        return POS_MAX
    return k


def parse_os_field(s: str):
    """Yield (metal, OS) pairs from a field like 'Cu26:3;Fe40:-1'."""
    s = (s or "").strip()
    if not s or s in ("N/A", "ERROR"):
        return
    for chunk in s.split(";"):
        m = ATOM_RE.match(chunk.strip())
        if not m:
            continue
        el = m.group(1)
        os_val = int(m.group(3))
        if el in TM_SET:
            yield el, os_val


def main():
    # First pass: collect all variants per (parent_dir, group_stem).
    groups = defaultdict(list)  # key -> [(charge, mult, basename, reactant, product)]

    n_rows = 0
    n_malformed = 0
    with open(IN_CSV, "r", newline="") as fh:
        r = csv.DictReader(fh)
        for row in r:
            n_rows += 1
            path = row["tar_zst_path"]
            parent, _, base = path.rpartition("/")
            parsed = parse_stem(base)
            if parsed is None:
                n_malformed += 1
                continue
            stem, charge, mult = parsed
            key = (parent, stem)
            groups[key].append((
                charge,
                mult,
                base,
                row["reactant_metal_oxidation_states"],
                row["product_metal_oxidation_states"],
            ))

    # Dedup: pick one variant per group following the charge/mult rule.
    picked = []
    n_skipped = 0
    charge_priority = {0: 0, -1: 1, 1: 2}

    for key, variants in groups.items():
        candidates = [v for v in variants if v[0] in charge_priority]
        if not candidates:
            n_skipped += 1
            continue
        candidates.sort(key=lambda v: (charge_priority[v[0]], v[1], v[2]))
        picked.append(candidates[0])

    # Aggregate reactant + product OS counts per metal.
    counts = defaultdict(Counter)  # metal -> Counter[os_bin]

    n_picked_with_data = 0
    for (_charge, _mult, _base, reactant, product) in picked:
        had_any = False
        for el, os_val in parse_os_field(reactant):
            counts[el][clamp(os_val)] += 1
            had_any = True
        for el, os_val in parse_os_field(product):
            counts[el][clamp(os_val)] += 1
            had_any = True
        if had_any:
            n_picked_with_data += 1

    # Write output CSV.
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        header = ["Metal"] + [str(k) for k in OS_COLS] + ["Total"]
        w.writerow(header)
        for metal in TM_ORDER:
            row_counts = counts.get(metal, Counter())
            row = [metal]
            for k in OS_COLS:
                row.append(row_counts.get(k, 0))
            row.append(sum(row_counts.values()))
            w.writerow(row)

    print(f"Input rows:                 {n_rows:,}", file=sys.stderr)
    print(f"Malformed paths:            {n_malformed:,}", file=sys.stderr)
    print(f"Reaction groups (pre-dedup):{len(groups):,}", file=sys.stderr)
    print(f"Groups skipped (no |c|<=1): {n_skipped:,}", file=sys.stderr)
    print(f"Picked variants:            {len(picked):,}", file=sys.stderr)
    print(f"  ...with TM atoms:         {n_picked_with_data:,}", file=sys.stderr)
    print(f"Wrote: {OUT_CSV}", file=sys.stderr)


if __name__ == "__main__":
    main()
