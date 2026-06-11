#!/usr/bin/env python3
"""Parallel version of extract_charge_spin.py.

Walks all xyz files under /groups/bsavoie2/zli43/Original-Files/doi_files/,
distributes line-2 parsing across N workers (multiprocessing.Pool), and
writes the same outputs as the serial script:

  results/charge_spin_groundtruth.csv
  results/charge_spin_coverage.csv
  results/charge_spin_value_distribution.csv
  results/charge_spin_unmatched_sample.txt
  results/charge_spin_summary.txt

Usage:
  python3 extract_charge_spin_parallel.py [N_WORKERS]
  default N_WORKERS = $(nproc) or 16, whichever is smaller.

Stage 1 (walk + collect xyz paths) is single-threaded.
Stage 2 (line-2 parse) is parallel across N workers.
Stage 3 (aggregate + write) is single-threaded again.
"""

import os
import re
import sys
import csv
from collections import Counter, defaultdict
from pathlib import Path
from multiprocessing import Pool, cpu_count
from time import time

DOI_ROOT = Path("/groups/bsavoie2/zli43/Original-Files/doi_files")
RESULTS  = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")


# ─────────────────────────────────────────────────────────────
# Path utilities
# ─────────────────────────────────────────────────────────────

def is_target_xyz(name):
    return name.endswith(".xyz") and not name.endswith("repacked.xyz")


def doi_from_path_parts(parts):
    for i, part in enumerate(parts):
        if part.startswith("10.") and i + 1 < len(parts):
            return part + "/" + parts[i + 1]
        if part == "no_doi":
            return "no_doi/" + "/".join(parts[i+1:i+3])
    return "?"


# ─────────────────────────────────────────────────────────────
# Pattern definitions — module-level so workers can import them
# ─────────────────────────────────────────────────────────────

_P1_CHARGE = re.compile(r"\b(?:charge|chg|q|total\s*charge)\s*[:=]?\s*([+-]?\d+)",
                        re.IGNORECASE)
_P1_MULT   = re.compile(r"\b(?:multiplicity|mult|2s\+1|mul)\s*[:=]?\s*(\d+)",
                        re.IGNORECASE)
_P1_MULT_M = re.compile(r"(?:^|[\s,;])M\s*[:=]\s*(\d+)\b")

_MULT_WORDS = {
    "singlet": 1, "doublet": 2, "triplet": 3, "quartet": 4,
    "quintet": 5, "sextet": 6, "septet": 7, "octet": 8,
}
_CHARGE_WORDS = {
    "neutral": 0,
    "anion": -1, "anionic": -1,
    "monoanion": -1, "monoanionic": -1,
    "cation": +1, "cationic": +1,
    "monocation": +1, "monocationic": +1,
    "dianion": -2, "dianionic": -2,
    "dication": +2, "dicationic": +2,
    "trianion": -3, "tricationic": +3, "trication": +3,
    "zwitterion": 0,
}

_P3        = re.compile(r"^\s*([+-]?\d+)\s+(\d+)\s*$")
_P3_LOOSE  = re.compile(r"^\s*([+-]?\d{1,2})\s+([1-9]\d?)\b")
_P4        = re.compile(r"\bS\s*=\s*(\d+(?:\.\d+)?|\d+/\d+)\b", re.IGNORECASE)


def parse_S_to_mult(s):
    if "/" in s:
        a, b = s.split("/")
        S = float(a) / float(b)
    else:
        S = float(s)
    mult = int(round(2*S + 1))
    return mult if 1 <= mult <= 10 else None


def parse_comment(line):
    """Return (charge, mult, pattern_label) or None."""
    if not line or not line.strip():
        return None
    L = line.strip()
    if len(L) > 500:
        L = L[:500]
    low = L.lower()

    c1 = _P1_CHARGE.search(L)
    m1 = _P1_MULT.search(L) or _P1_MULT_M.search(L)
    if c1 and m1:
        try:
            ch = int(c1.group(1))
            mu = int(m1.group(1))
            if -5 <= ch <= 5 and 1 <= mu <= 10:
                return (ch, mu, "P1_explicit")
        except ValueError:
            pass
    p1_ch = None; p1_mu = None
    if c1:
        try:
            v = int(c1.group(1))
            if -5 <= v <= 5: p1_ch = v
        except ValueError:
            pass
    if m1:
        try:
            v = int(m1.group(1))
            if 1 <= v <= 10: p1_mu = v
        except ValueError:
            pass

    p2_ch = None; p2_mu = None
    for w, val in _MULT_WORDS.items():
        if re.search(r"\b" + w + r"\b", low):
            p2_mu = val
            break
    for w, val in _CHARGE_WORDS.items():
        if re.search(r"\b" + w + r"\b", low):
            p2_ch = val
            break

    ch_candidate = p1_ch if p1_ch is not None else p2_ch
    mu_candidate = p1_mu if p1_mu is not None else p2_mu

    if ch_candidate is not None and mu_candidate is not None:
        label = ("P1_explicit" if (p1_ch is not None and p1_mu is not None)
                 else "P1+P2"  if (p1_ch is not None or p1_mu is not None)
                 else "P2_words")
        return (ch_candidate, mu_candidate, label)

    m3 = _P3.match(L)
    if m3:
        ch = int(m3.group(1)); mu = int(m3.group(2))
        if -5 <= ch <= 5 and 1 <= mu <= 10:
            return (ch, mu, "P3_two_ints_strict")

    m4 = _P4.search(L)
    if m4:
        mu_p4 = parse_S_to_mult(m4.group(1))
        if mu_p4 is not None:
            ch_use = ch_candidate if ch_candidate is not None else 0
            label = ("P4_spin+P1" if p1_ch is not None
                     else "P4_spin+P2" if p2_ch is not None
                     else "P4_spin_only")
            return (ch_use, mu_p4, label)

    m3l = _P3_LOOSE.match(L)
    if m3l:
        ch = int(m3l.group(1)); mu = int(m3l.group(2))
        if -5 <= ch <= 5 and 1 <= mu <= 10:
            if len(L) < 80 and not any(w in low for w in ("=", "//", "transition", "trans-state", "atom")):
                return (ch, mu, "P3_two_ints_loose")

    if ch_candidate is not None and mu_candidate is None:
        return (ch_candidate, None, "partial_charge_only")
    if mu_candidate is not None and ch_candidate is None:
        return (None, mu_candidate, "partial_mult_only")

    return None


# ─────────────────────────────────────────────────────────────
# Worker function (runs in each child process)
# ─────────────────────────────────────────────────────────────

def process_one(path_str):
    """Read line 2 of path; return tuple (path_str, doi, leaf_dir,
       charge, mult, pattern, comment_excerpt, status)
       where status is one of: 'blank', 'matched', 'partial', 'unmatched', 'unreadable'.
       For 'unmatched' we still return the comment text so the master can sample it."""
    p = Path(path_str)
    parts = p.parts
    doi = doi_from_path_parts(parts)
    leaf = str(p.parent)
    try:
        with open(path_str, errors="replace") as f:
            f.readline()
            comment = f.readline()
    except OSError:
        return (path_str, doi, leaf, None, None, None, "", "unreadable")

    if comment is None or not comment.strip():
        return (path_str, doi, leaf, None, None, None, "", "blank")

    comment_stripped = comment.rstrip("\n\r").strip()
    excerpt = comment_stripped[:200]

    parsed = parse_comment(comment_stripped)
    if parsed is None:
        return (path_str, doi, leaf, None, None, None, excerpt, "unmatched")

    ch, mu, label = parsed
    status = "partial" if label.startswith("partial") else "matched"
    return (path_str, doi, leaf, ch, mu, label, excerpt, status)


# ─────────────────────────────────────────────────────────────
# Master / aggregation
# ─────────────────────────────────────────────────────────────

def collect_xyz_paths():
    """Walk doi_files; return list of leaf xyz file path strings."""
    paths = []
    for dp, dnames, fnames in os.walk(DOI_ROOT):
        if dnames:
            continue
        for fn in fnames:
            if is_target_xyz(fn):
                paths.append(os.path.join(dp, fn))
    return paths


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else min(cpu_count(), 16)

    t0 = time()
    print(f"[{0:.0f}s] Stage 1: walking filesystem ...", file=sys.stderr)
    paths = collect_xyz_paths()
    print(f"[{time()-t0:.0f}s]   {len(paths):,} xyz files queued for parsing",
          file=sys.stderr)

    print(f"[{time()-t0:.0f}s] Stage 2: parsing with {n_workers} workers ...",
          file=sys.stderr)

    # Aggregators
    pattern_counter = Counter()
    charge_counter  = Counter()
    mult_counter    = Counter()
    # per-paper stats: doi -> {"matched": int, "total": int,
    #                          "all_leaves": set, "leaves_with_label": set}
    pp = defaultdict(lambda: {"matched": 0, "total": 0,
                              "all_leaves": set(),
                              "leaves_with_label": set()})
    unmatched_sample = []
    UNMATCHED_LIMIT = 500

    n_total = 0; n_blank = 0; n_matched = 0; n_partial = 0; n_unreadable = 0

    gt_path = RESULTS / "charge_spin_groundtruth.csv"
    gt_fp = open(gt_path, "w", newline="")
    gt_w = csv.writer(gt_fp)
    gt_w.writerow(["file", "doi", "leaf_dir", "charge", "multiplicity",
                   "pattern", "raw_comment"])

    last_print = 0
    with Pool(n_workers) as pool:
        # chunk size: bigger = less overhead, but worse load balance
        for res in pool.imap_unordered(process_one, paths, chunksize=400):
            path, doi, leaf, ch, mu, label, excerpt, status = res
            n_total += 1
            pp[doi]["total"] += 1
            pp[doi]["all_leaves"].add(leaf)

            if status == "blank":
                n_blank += 1
            elif status == "unreadable":
                n_unreadable += 1
            elif status == "unmatched":
                if len(unmatched_sample) < UNMATCHED_LIMIT:
                    unmatched_sample.append((path, excerpt))
            else:  # matched or partial
                pattern_counter[label] += 1
                if ch is not None:
                    charge_counter[ch] += 1
                if mu is not None:
                    mult_counter[mu] += 1
                if status == "partial":
                    n_partial += 1
                else:
                    n_matched += 1
                    pp[doi]["matched"] += 1
                    pp[doi]["leaves_with_label"].add(leaf)
                gt_w.writerow([path, doi, leaf,
                               "" if ch is None else ch,
                               "" if mu is None else mu,
                               label, excerpt])

            if n_total - last_print >= 100000:
                print(f"[{time()-t0:.0f}s]   {n_total:,}/{len(paths):,} "
                      f"matched={n_matched:,}  partial={n_partial:,}  "
                      f"unread={n_unreadable:,}", file=sys.stderr)
                last_print = n_total

    gt_fp.close()
    print(f"[{time()-t0:.0f}s] Stage 2 done. Writing summary files ...",
          file=sys.stderr)

    # ── coverage CSV
    cov_path = RESULTS / "charge_spin_coverage.csv"
    with open(cov_path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["doi", "n_xyz_total", "n_xyz_matched",
                    "n_leaves_total", "n_leaves_with_any_label"])
        for doi in sorted(pp):
            s = pp[doi]
            w.writerow([doi, s["total"], s["matched"],
                        len(s["all_leaves"]),
                        len(s["leaves_with_label"])])

    # ── value distribution
    val_path = RESULTS / "charge_spin_value_distribution.csv"
    with open(val_path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["dimension", "value", "count"])
        for ch in sorted(charge_counter):
            w.writerow(["charge", ch, charge_counter[ch]])
        for mu in sorted(mult_counter):
            w.writerow(["multiplicity", mu, mult_counter[mu]])
        for pat, c in pattern_counter.most_common():
            w.writerow(["pattern", pat, c])

    # ── unmatched sample
    um_path = RESULTS / "charge_spin_unmatched_sample.txt"
    with open(um_path, "w") as fp:
        fp.write(f"First {len(unmatched_sample)} non-blank comment lines that "
                 f"did NOT match any pattern.\n\n")
        for path, line in unmatched_sample:
            fp.write(f"{path}\n  {line!r}\n\n")

    # ── summary
    n_papers_total       = len(pp)
    n_papers_with_label  = sum(1 for d in pp.values() if d["matched"] > 0)
    n_leaves_total       = sum(len(d["all_leaves"]) for d in pp.values())
    n_leaves_with_label  = sum(len(d["leaves_with_label"]) for d in pp.values())

    sum_path = RESULTS / "charge_spin_summary.txt"
    with open(sum_path, "w") as fp:
        fp.write("Charge / multiplicity extraction summary (parallel run)\n")
        fp.write("=" * 60 + "\n\n")
        fp.write(f"workers used: {n_workers}\n")
        fp.write(f"elapsed wallclock: {time()-t0:.1f} seconds\n\n")
        fp.write(f"Total xyz files scanned:                 {n_total:,}\n")
        fp.write(f"  blank or unreadable comment lines:     {n_blank+n_unreadable:,}\n")
        fp.write(f"    of which unreadable:                 {n_unreadable:,}\n")
        fp.write(f"  matched (full charge+mult):            {n_matched:,}  "
                 f"({100*n_matched/max(n_total,1):.2f}%)\n")
        fp.write(f"  partial (charge only OR mult only):    {n_partial:,}  "
                 f"({100*n_partial/max(n_total,1):.2f}%)\n")
        fp.write(f"  unmatched (non-blank but no signal):   "
                 f"{n_total - n_matched - n_partial - n_blank - n_unreadable:,}\n\n")

        fp.write(f"Pattern frequency:\n")
        for pat, c in pattern_counter.most_common():
            fp.write(f"  {pat:30s} {c:>9,}\n")

        fp.write(f"\n--- per-paper coverage ---\n")
        fp.write(f"Total papers seen:                   {n_papers_total:,}\n")
        fp.write(f"Papers with at least 1 labeled xyz:  {n_papers_with_label:,}  "
                 f"({100*n_papers_with_label/max(n_papers_total,1):.2f}%)\n")
        fp.write(f"Total leaf dirs:                     {n_leaves_total:,}\n")
        fp.write(f"Leaves with at least 1 labeled xyz:  {n_leaves_with_label:,}  "
                 f"({100*n_leaves_with_label/max(n_leaves_total,1):.2f}%)\n\n")

        fp.write(f"--- charge distribution ---\n")
        for ch in sorted(charge_counter):
            fp.write(f"  charge={ch:+3d}    {charge_counter[ch]:>9,}\n")
        fp.write(f"\n--- multiplicity distribution ---\n")
        for mu in sorted(mult_counter):
            fp.write(f"  M={mu}    {mult_counter[mu]:>9,}\n")
        fp.write(f"\nOutputs:\n  {gt_path}\n  {cov_path}\n  {val_path}\n  "
                 f"{um_path}\n  {sum_path}\n")

    print(f"[{time()-t0:.0f}s] Done. Outputs in {RESULTS}", file=sys.stderr)


if __name__ == "__main__":
    main()
