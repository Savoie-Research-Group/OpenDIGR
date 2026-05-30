#!/usr/bin/env python3
"""
audit_metal_os.py

Per-element oxidation-state audit for the radial-dial panel.

Generalization of audit_iridium_os.py. The dial bins clamp OS >= +6 into the
+6 column (see build_tm_os_matrix.py), so the bar may bundle genuine +6 atoms
together with extractor-artifact values at +7, +8, +9. This script re-runs
the v2 dedup on transition_metal_oxidation_states.csv for one or more target
elements and writes:

  <el>_os_summary.txt         per-OS atom counts (reactant / product / total),
                              unclamped and clamped (the form the dial uses).
  <el>_reactions.csv          one row per kept reaction containing the element,
                              with atoms listed as
                              "<el><idx>:<reactant_os>->,<product_os>".
  <el>_high_os_reactions.csv  subset where any target atom has OS >= +6.

Usage:
  python3 audit_metal_os.py Ir Pt W Re Os
  python3 audit_metal_os.py --all   # Ir Pt W Re Os
"""

import csv
import os
import re
import sys
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
IN_CSV = os.path.join(ROOT, "transition_metal_oxidation_states.csv")

DEFAULT_ELEMENTS = ["Ir", "Pt", "W", "Re", "Os"]
HIGH_OS_THRESHOLD = 6  # matches POS_MAX clamp in build_tm_os_matrix.py
CLAMP_MIN, CLAMP_MAX = -3, 6

ATOM_RE = re.compile(r"([A-Z][a-z]?)(\d+):(-?\d+)")


def parse_stem(basename):
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


def parse_atoms(field):
    field = (field or "").strip()
    if not field or field in ("N/A", "ERROR"):
        return []
    out = []
    for chunk in field.split(";"):
        m = ATOM_RE.match(chunk.strip())
        if m:
            out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return out


def clamp(k):
    if k <= CLAMP_MIN:
        return CLAMP_MIN
    if k >= CLAMP_MAX:
        return CLAMP_MAX
    return k


def load_picked():
    """Return (n_rows, picked_variants) using the v2 dedup rule."""
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
    return n_rows, len(groups), picked


def audit_element(target_el, n_rows, n_groups, picked):
    reactant_os = Counter()
    product_os = Counter()
    kept = []

    for charge, mult, base, path, r_field, p_field in picked:
        r_atoms = [a for a in parse_atoms(r_field) if a[0] == target_el]
        p_atoms = [a for a in parse_atoms(p_field) if a[0] == target_el]
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
            f"{target_el}{i}:{r_by_idx.get(i, '?')}->{p_by_idx.get(i, '?')}"
            for i in all_idx
        )
        max_os = max(
            [os for _, _, os in r_atoms] + [os for _, _, os in p_atoms]
        )
        kept.append({
            "path": path, "charge": charge, "mult": mult,
            "n_atoms": len(all_idx), "max_os": max_os, "atoms": atoms_repr,
        })

    all_os = sorted(set(reactant_os) | set(product_os))
    raw_total = Counter({k: reactant_os[k] + product_os[k] for k in all_os})
    clamped_total = Counter()
    for k, n in raw_total.items():
        clamped_total[clamp(k)] += n

    out_summary = os.path.join(HERE, f"{target_el.lower()}_os_summary.txt")
    out_reactions = os.path.join(HERE, f"{target_el.lower()}_reactions.csv")
    out_high = os.path.join(HERE, f"{target_el.lower()}_high_os_reactions.csv")

    with open(out_summary, "w") as fh:
        fh.write(f"{target_el} oxidation-state audit\n")
        fh.write(f"Source CSV : {IN_CSV}\n")
        fh.write(f"Rows read  : {n_rows:,}\n")
        fh.write(f"Reaction groups (pre-dedup): {n_groups:,}\n")
        fh.write(f"Picked variants (post-dedup): {len(picked):,}\n")
        fh.write(f"Reactions containing {target_el}: {len(kept):,}\n\n")

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
        if sane_total:
            pct = 100 * high_atoms / sane_total
        else:
            pct = 0.0
        fh.write(
            f"\nAtoms with OS >= +{HIGH_OS_THRESHOLD}: "
            f"{high_atoms:,} / {sane_total:,} ({pct:.1f}%) "
            f"-- all collapsed into the +6 dial bin.\n"
        )

    fields = ["tar_zst_path", "charge", "multiplicity",
              f"n_{target_el.lower()}_atoms",
              f"max_{target_el.lower()}_os",
              f"{target_el.lower()}_atoms_reactant->product"]
    with open(out_reactions, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in sorted(kept, key=lambda r: (-r["max_os"], r["path"])):
            w.writerow([r["path"], r["charge"], r["mult"],
                        r["n_atoms"], r["max_os"], r["atoms"]])

    high_rows = [r for r in kept if r["max_os"] >= HIGH_OS_THRESHOLD]
    with open(out_high, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(fields)
        for r in sorted(high_rows, key=lambda r: (-r["max_os"], r["path"])):
            w.writerow([r["path"], r["charge"], r["mult"],
                        r["n_atoms"], r["max_os"], r["atoms"]])

    return {
        "el": target_el, "reactions": len(kept),
        "atoms_total": sane_total, "atoms_high": high_atoms,
        "pct_high": pct, "raw_total": raw_total,
        "summary": out_summary, "reactions_csv": out_reactions, "high_csv": out_high,
    }


def main(argv):
    if len(argv) > 1 and argv[1] == "--all":
        elements = DEFAULT_ELEMENTS
    elif len(argv) > 1:
        elements = argv[1:]
    else:
        elements = DEFAULT_ELEMENTS

    print("Loading + deduping CSV ...", file=sys.stderr)
    n_rows, n_groups, picked = load_picked()
    print(f"  rows={n_rows:,}  groups={n_groups:,}  picked={len(picked):,}",
          file=sys.stderr)

    print(f"\n{'el':>3} | {'reactions':>10} | {'atoms':>9} | "
          f"{'>=+6 atoms':>10} | {'%':>6} | top OS atoms (raw, top 5)")
    print("-" * 96)
    for el in elements:
        info = audit_element(el, n_rows, n_groups, picked)
        top = sorted(info["raw_total"].items(),
                     key=lambda kv: -kv[1])[:5]
        top_str = " ".join(f"{k:+d}:{n}" for k, n in top)
        print(f"{el:>3} | {info['reactions']:>10,} | {info['atoms_total']:>9,} | "
              f"{info['atoms_high']:>10,} | {info['pct_high']:>5.1f}% | {top_str}")
        print(f"    -> {info['summary']}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv)
