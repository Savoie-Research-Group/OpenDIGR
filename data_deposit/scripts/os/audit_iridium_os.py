#!/usr/bin/env python3
"""
audit_iridium_os.py

Iridium oxidation-state audit for the radial-dial panel.

A reviewer flagged the unusually large +6 bin for late TMs (Ir, Pt, W, Re, Os).
The dial bins clamp every OS >= +6 into the +6 column (see build_tm_os_matrix.py),
so the bar bundles together genuine +6 atoms and any high-OS values the
extractor returned (+7, +8, +9, ...). This script re-runs the v2 dedup on
transition_metal_oxidation_states.csv, isolates iridium, and writes:

  ir_os_summary.txt          per-OS atom counts (reactant / product / total),
                             both unclamped and the clamped form that matches
                             the radial-dial bin.
  ir_reactions.csv           one row per kept reaction that contains Ir,
                             with the Ir atoms listed as
                             "<index>:<reactant_os>->,<product_os>".
  ir_high_os_reactions.csv   subset where any Ir atom has OS >= 6.

Dedup matches build_tm_os_matrix.py: per (parent_dir, stem-minus-last-two-_-tokens),
prefer charge 0, else |charge|==1 (-1 preferred over +1), then lowest multiplicity,
then alphabetical basename.
"""

import csv
import os
import re
import sys
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
IN_CSV = os.path.join(ROOT, "transition_metal_oxidation_states.csv")

OUT_SUMMARY = os.path.join(HERE, "ir_os_summary.txt")
OUT_REACTIONS = os.path.join(HERE, "ir_reactions.csv")
OUT_HIGH = os.path.join(HERE, "ir_high_os_reactions.csv")

TARGET_EL = "Ir"
HIGH_OS_THRESHOLD = 6  # matches POS_MAX clamp in build_tm_os_matrix.py
CLAMP_MIN, CLAMP_MAX = -3, 6

ATOM_RE = re.compile(r"([A-Z][a-z]?)(\d+):(-?\d+)")


def parse_stem(basename: str):
    if not basename.endswith(".tar.zst"):
        return None
    stem = basename[:-len(".tar.zst")]
    parts = stem.rsplit("_", 2)
    if len(parts) != 3:
        return None
    try:
        return parts[0], int(parts[1]), int(parts[2])
    except ValueError:
        return None


def parse_atoms(field: str):
    field = (field or "").strip()
    if not field or field in ("N/A", "ERROR"):
        return []
    out = []
    for chunk in field.split(";"):
        m = ATOM_RE.match(chunk.strip())
        if m:
            out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return out


def clamp(k: int) -> int:
    if k <= CLAMP_MIN:
        return CLAMP_MIN
    if k >= CLAMP_MAX:
        return CLAMP_MAX
    return k


def main():
    groups = defaultdict(list)
    n_rows = 0
    with open(IN_CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            n_rows += 1
            path = row["tar_zst_path"]
            parent, _, base = path.rpartition("/")
            parsed = parse_stem(base)
            if parsed is None:
                continue
            stem, charge, mult = parsed
            groups[(parent, stem)].append((
                charge, mult, base, path,
                row["reactant_metal_oxidation_states"],
                row["product_metal_oxidation_states"],
            ))

    charge_priority = {0: 0, -1: 1, 1: 2}
    picked = []
    for variants in groups.values():
        cands = [v for v in variants if v[0] in charge_priority]
        if not cands:
            continue
        cands.sort(key=lambda v: (charge_priority[v[0]], v[1], v[2]))
        picked.append(cands[0])

    reactant_os = Counter()
    product_os = Counter()
    ir_kept = []

    for charge, mult, base, path, r_field, p_field in picked:
        r_atoms = [a for a in parse_atoms(r_field) if a[0] == TARGET_EL]
        p_atoms = [a for a in parse_atoms(p_field) if a[0] == TARGET_EL]
        if not r_atoms and not p_atoms:
            continue
        for _, _, os_val in r_atoms:
            reactant_os[os_val] += 1
        for _, _, os_val in p_atoms:
            product_os[os_val] += 1
        r_by_idx = {idx: os for _, idx, os in r_atoms}
        p_by_idx = {idx: os for _, idx, os in p_atoms}
        all_idx = sorted(set(r_by_idx) | set(p_by_idx))
        atoms_repr = ";".join(
            f"Ir{i}:{r_by_idx.get(i, '?')}->{p_by_idx.get(i, '?')}"
            for i in all_idx
        )
        max_os = max(
            [os for _, _, os in r_atoms] + [os for _, _, os in p_atoms]
        )
        ir_kept.append({
            "path": path,
            "charge": charge,
            "mult": mult,
            "n_ir_atoms": len(all_idx),
            "max_ir_os": max_os,
            "ir_atoms": atoms_repr,
        })

    # --- ir_os_summary.txt -------------------------------------------------
    all_os = sorted(set(reactant_os) | set(product_os))
    raw_total = Counter()
    for k in all_os:
        raw_total[k] = reactant_os[k] + product_os[k]

    clamped_total = Counter()
    for k, n in raw_total.items():
        clamped_total[clamp(k)] += n

    with open(OUT_SUMMARY, "w") as fh:
        fh.write(f"Iridium oxidation-state audit\n")
        fh.write(f"Source CSV : {IN_CSV}\n")
        fh.write(f"Rows read  : {n_rows:,}\n")
        fh.write(f"Reaction groups (pre-dedup): {len(groups):,}\n")
        fh.write(f"Picked variants (post-dedup): {len(picked):,}\n")
        fh.write(f"Reactions containing Ir    : {len(ir_kept):,}\n\n")

        fh.write("Per-OS atom counts (deduped, unclamped):\n")
        fh.write(f"{'OS':>4}  {'reactant':>10}  {'product':>10}  {'total':>10}\n")
        for k in all_os:
            fh.write(f"{k:>+4d}  {reactant_os[k]:>10,}  "
                     f"{product_os[k]:>10,}  {raw_total[k]:>10,}\n")
        fh.write(f"{'sum':>4}  {sum(reactant_os.values()):>10,}  "
                 f"{sum(product_os.values()):>10,}  {sum(raw_total.values()):>10,}\n")

        fh.write("\nClamped to dial bins [-3..+6] (matches panel_b_dials.svg):\n")
        fh.write(f"{'bin':>4}  {'total':>10}\n")
        for k in range(CLAMP_MIN, CLAMP_MAX + 1):
            fh.write(f"{k:>+4d}  {clamped_total[k]:>10,}\n")

        high_atoms = sum(n for k, n in raw_total.items() if k >= HIGH_OS_THRESHOLD)
        sane_total = sum(raw_total.values())
        fh.write(
            f"\nAtoms with OS >= +{HIGH_OS_THRESHOLD}: "
            f"{high_atoms:,} / {sane_total:,} "
            f"({100*high_atoms/sane_total:.1f}%) -- all collapsed into the +6 dial bin.\n"
        )

    # --- ir_reactions.csv --------------------------------------------------
    fields = ["tar_zst_path", "charge", "multiplicity",
              "n_ir_atoms", "max_ir_os", "ir_atoms_reactant->product"]
    with open(OUT_REACTIONS, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in sorted(ir_kept, key=lambda r: (-r["max_ir_os"], r["path"])):
            w.writerow([r["path"], r["charge"], r["mult"],
                        r["n_ir_atoms"], r["max_ir_os"], r["ir_atoms"]])

    # --- ir_high_os_reactions.csv -----------------------------------------
    high_rows = [r for r in ir_kept if r["max_ir_os"] >= HIGH_OS_THRESHOLD]
    with open(OUT_HIGH, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in sorted(high_rows, key=lambda r: (-r["max_ir_os"], r["path"])):
            w.writerow([r["path"], r["charge"], r["mult"],
                        r["n_ir_atoms"], r["max_ir_os"], r["ir_atoms"]])

    print(f"Rows read         : {n_rows:,}", file=sys.stderr)
    print(f"Reaction groups   : {len(groups):,}", file=sys.stderr)
    print(f"Picked variants   : {len(picked):,}", file=sys.stderr)
    print(f"Ir-bearing kept   : {len(ir_kept):,}", file=sys.stderr)
    print(f"  ... max OS >= +{HIGH_OS_THRESHOLD}: {len(high_rows):,}",
          file=sys.stderr)
    print(f"Wrote: {OUT_SUMMARY}", file=sys.stderr)
    print(f"Wrote: {OUT_REACTIONS}", file=sys.stderr)
    print(f"Wrote: {OUT_HIGH}", file=sys.stderr)


if __name__ == "__main__":
    main()
