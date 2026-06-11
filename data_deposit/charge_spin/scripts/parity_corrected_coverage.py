#!/usr/bin/env python3
"""Recompute lowest-spin coverage with parity correction.

For each full-match row in charge_spin_groundtruth.csv (charge AND
multiplicity both labeled), this script:

  1. Opens the xyz file and parses element symbols from the
     coordinate lines.
  2. Computes total electron count = sum(Z) - charge.
  3. Determines the system-specific "lowest possible" multiplicity:
        even electron count → M_min = 1   (closed-shell singlet)
        odd electron count  → M_min = 2   (one unpaired electron)
  4. Compares the system-specific M_min to the author-labeled M.
  5. Computes coverage of progressive sampling depths
     {M_min}, {M_min, M_min+2}, {M_min, M_min+2, M_min+4}, ...
     (because spin states change in steps of 2 within a fixed
     electron-count parity).

Writes:
  results/parity_corrected_coverage.csv     row per file with
                                            atoms, Z_sum, electron
                                            count, parity, M_min,
                                            labeled M, match flag
  results/parity_corrected_summary.txt      coverage tables
"""

import csv
import re
import sys
from collections import Counter
from pathlib import Path

GT_CSV   = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results/charge_spin_groundtruth.csv")
OUT_CSV  = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results/parity_corrected_coverage.csv")
OUT_TXT  = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results/parity_corrected_summary.txt")

# Atomic numbers Z = 1..118 by symbol
_PT = """
H  He
Li Be B  C  N  O  F  Ne
Na Mg Al Si P  S  Cl Ar
K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr
Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
""".split()
Z_BY_SYMBOL = {sym: i + 1 for i, sym in enumerate(_PT)}

# Same coordinate regex as element extractor
_NUM = r"-?\d+(\.\d+)?([eE][-+]?\d+)?"
_COORD_RE = re.compile(r"^\s*([A-Z][a-z]?)\s+({n})\s+({n})\s+({n})".format(n=_NUM))


def parse_atoms_z_sum(path, max_lines=20000):
    """Walk xyz file; sum Z over all atom lines we can identify."""
    z_sum = 0
    n_atoms = 0
    n_unknown = 0
    try:
        with open(path, errors="replace") as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    break
                m = _COORD_RE.match(line)
                if not m:
                    continue
                sym = m.group(1)
                z = Z_BY_SYMBOL.get(sym)
                if z is None:
                    n_unknown += 1
                    continue
                z_sum += z
                n_atoms += 1
    except OSError:
        return (0, 0, 0)
    return (n_atoms, z_sum, n_unknown)


def cumulative_coverage(actual_M, M_min):
    """Return the smallest depth k such that {M_min, M_min+2, ..., M_min+2k}
    contains actual_M.  None if actual_M is in the wrong parity track or
    below M_min."""
    if actual_M < M_min:
        return None
    if (actual_M - M_min) % 2 != 0:
        # mismatched parity: author's M is wrong-parity for our Z-derived count
        return None
    return (actual_M - M_min) // 2


def main():
    print("Reading ground-truth CSV ...", file=sys.stderr)

    n_rows = 0
    n_full_match = 0
    n_processed = 0
    n_unreadable = 0
    n_parity_mismatch = 0

    coverage_depth = Counter()   # depth k -> count
    parity_dist    = Counter()   # "even"/"odd"
    mmin_dist      = Counter()
    deltaM_dist    = Counter()   # actual_M - M_min (signed)

    # also: per-pattern accuracy (do high-confidence tiers agree more?)
    by_pattern = {}

    with open(OUT_CSV, "w", newline="") as out_fp:
        wout = csv.writer(out_fp)
        wout.writerow(["file", "doi", "leaf_dir", "charge", "labeled_M",
                       "pattern", "n_atoms", "z_sum", "electron_count",
                       "parity", "M_min", "depth_k", "matched"])

        with open(GT_CSV) as fp:
            r = csv.DictReader(fp)
            for row in r:
                n_rows += 1
                pat = row["pattern"]
                # need both charge and multiplicity to do parity check
                if pat.startswith("partial"):
                    continue
                if row["charge"] == "" or row["multiplicity"] == "":
                    continue
                try:
                    ch = int(row["charge"])
                    mu = int(row["multiplicity"])
                except ValueError:
                    continue
                n_full_match += 1
                if n_full_match % 10000 == 0:
                    print(f"  ... {n_full_match} full-match rows processed",
                          file=sys.stderr)

                n_atoms, z_sum, _ = parse_atoms_z_sum(row["file"])
                if n_atoms == 0:
                    n_unreadable += 1
                    continue
                n_processed += 1

                electrons = z_sum - ch
                parity = "even" if electrons % 2 == 0 else "odd"
                M_min = 1 if parity == "even" else 2

                parity_dist[parity] += 1
                mmin_dist[M_min] += 1
                deltaM_dist[mu - M_min] += 1

                depth = cumulative_coverage(mu, M_min)
                matched = depth is not None
                if depth is not None:
                    coverage_depth[depth] += 1
                else:
                    n_parity_mismatch += 1

                # per-pattern
                key = pat
                if key not in by_pattern:
                    by_pattern[key] = {"total": 0, "matched_at_k0": 0,
                                       "matched_at_k1": 0, "matched_at_k2": 0,
                                       "parity_mismatch": 0}
                bp = by_pattern[key]
                bp["total"] += 1
                if depth is None:
                    bp["parity_mismatch"] += 1
                else:
                    if depth <= 0: bp["matched_at_k0"] += 1
                    if depth <= 1: bp["matched_at_k1"] += 1
                    if depth <= 2: bp["matched_at_k2"] += 1

                wout.writerow([row["file"], row["doi"], row["leaf_dir"],
                               ch, mu, pat,
                               n_atoms, z_sum, electrons,
                               parity, M_min,
                               "" if depth is None else depth,
                               "1" if matched else "0"])

    print(f"Done. {n_processed} processed, {n_unreadable} unreadable, "
          f"{n_parity_mismatch} parity mismatches.", file=sys.stderr)

    # ── compute cumulative coverage curves
    total_for_pct = n_processed
    cum = 0
    cum_rows = []
    for k in sorted(coverage_depth):
        cum += coverage_depth[k]
        cum_rows.append((k, coverage_depth[k], cum, 100*cum/total_for_pct))

    # ── write summary
    with open(OUT_TXT, "w") as fp:
        fp.write("Parity-corrected lowest-spin coverage summary\n")
        fp.write("=" * 60 + "\n\n")
        fp.write(f"Total ground-truth rows (full matches): {n_full_match:,}\n")
        fp.write(f"  successfully parsed (xyz readable):    {n_processed:,}\n")
        fp.write(f"  unreadable / empty xyz coordinates:    {n_unreadable:,}\n")
        fp.write(f"  parity mismatch (author M wrong parity for our Z-sum): {n_parity_mismatch:,}  "
                 f"({100*n_parity_mismatch/max(n_processed,1):.2f}%)\n\n")

        fp.write("--- Electron-count parity distribution ---\n")
        for p in ("even", "odd"):
            fp.write(f"  {p:5s}    {parity_dist[p]:>9,}  "
                     f"({100*parity_dist[p]/max(n_processed,1):.1f}%)\n")
        fp.write("\n--- System-specific M_min distribution ---\n")
        for m in sorted(mmin_dist):
            fp.write(f"  M_min={m}   {mmin_dist[m]:>9,}  "
                     f"({100*mmin_dist[m]/max(n_processed,1):.1f}%)\n")

        fp.write("\n--- Author M minus system-specific M_min ---\n")
        fp.write("(positive = author chose a higher spin state than the minimum)\n")
        fp.write("(odd values indicate a parity mismatch between author M and Z-derived parity)\n")
        for d in sorted(deltaM_dist):
            tag = "  (parity OK)" if d % 2 == 0 else "  ** parity mismatch **"
            fp.write(f"  ΔM={d:+d}    {deltaM_dist[d]:>9,}{tag}\n")

        fp.write("\n--- Cumulative coverage of 'sample M_min, M_min+2, ..., M_min+2k' ---\n")
        fp.write(f"{'depth k':>8}  {'M tried':>14}  {'count':>9}  {'cumulative':>10}  {'cum %':>7}\n")
        for k, c, cum_c, pct in cum_rows:
            mlist = "M_min" if k == 0 else f"+{k*2} from M_min"
            fp.write(f"{k:>8}  {mlist:>14}  {c:>9,}  {cum_c:>10,}  {pct:>6.2f}%\n")
        # note: parity_mismatch rows are excluded from cum_rows by construction
        fp.write(f"\nNote: {n_parity_mismatch:,} parity mismatches "
                 f"({100*n_parity_mismatch/max(n_processed,1):.2f}%) are not included "
                 "in the cumulative coverage above — they cannot be reached at any\n"
                 "spin-state sampling depth that respects electron-count parity.\n")

        fp.write("\n--- Per-pattern coverage ---\n")
        fp.write(f"{'pattern':<32s} {'N':>9s} {'M_min':>8s} {'+2':>8s} {'+4':>8s} {'mismatch':>10s}\n")
        for key in sorted(by_pattern, key=lambda k: -by_pattern[k]["total"]):
            bp = by_pattern[key]
            t = bp["total"]
            if t == 0: continue
            fp.write(f"{key:<32s} {t:>9,} "
                     f"{100*bp['matched_at_k0']/t:>7.1f}% "
                     f"{100*bp['matched_at_k1']/t:>7.1f}% "
                     f"{100*bp['matched_at_k2']/t:>7.1f}% "
                     f"{100*bp['parity_mismatch']/t:>9.1f}%\n")

        fp.write(f"\nOutputs:\n  {OUT_CSV}\n  {OUT_TXT}\n")


if __name__ == "__main__":
    main()
