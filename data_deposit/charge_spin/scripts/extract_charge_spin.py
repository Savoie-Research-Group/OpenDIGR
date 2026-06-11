#!/usr/bin/env python3
"""Extract explicit charge/multiplicity labels from xyz comment lines.

Walks every *.xyz under /groups/bsavoie2/zli43/Original-Files/doi_files/,
reads only line 2 (the standard xyz comment line), and applies a stack
of regex patterns to detect explicit charge / multiplicity / word forms.

Outputs (all in results/):
  charge_spin_groundtruth.csv      one row per xyz file that matched
  charge_spin_coverage.csv         per-paper / per-leaf coverage summary
  charge_spin_value_distribution.csv  histogram of (charge, mult)
  charge_spin_unmatched_sample.txt  first 500 non-empty unmatched comments
  charge_spin_summary.txt          human-readable summary

Patterns, in order of confidence (a higher one wins; lower ones are
only tried if higher ones fail):

  P1  "charge=...  mult=..."        explicit named fields
  P2  word forms (singlet/doublet/triplet, neutral/anion/cation, …)
  P3  "two ints" on the comment line ("0 1", "-1 2") — Gaussian header style
  P4  "S = 0", "S = 1/2"            spin quantum number; mult = 2S+1
"""

import os
import re
import sys
import csv
from collections import Counter, defaultdict
from pathlib import Path

DOI_ROOT = Path("/groups/bsavoie2/zli43/Original-Files/doi_files")
RESULTS  = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")

# Skip *repacked.xyz to be consistent with prior counts.
def is_target_xyz(name):
    return name.endswith(".xyz") and not name.endswith("repacked.xyz")


def doi_from_path(p):
    parts = p.parts
    for i, part in enumerate(parts):
        if part.startswith("10.") and i + 1 < len(parts):
            return part + "/" + parts[i + 1]
        if part == "no_doi":
            return "no_doi/" + "/".join(parts[i+1:i+3])
    return None


# ─────────────────────────────────────────────────────────────
# Pattern definitions
# ─────────────────────────────────────────────────────────────

# P1 — explicit named fields anywhere in the line
_P1_CHARGE = re.compile(r"\b(?:charge|chg|q|total\s*charge)\s*[:=]?\s*([+-]?\d+)",
                        re.IGNORECASE)
_P1_MULT   = re.compile(r"\b(?:multiplicity|mult|2s\+1|mul)\s*[:=]?\s*(\d+)",
                        re.IGNORECASE)
# also allow "M=2" but require it preceded by space/start or after non-alpha
_P1_MULT_M = re.compile(r"(?:^|[\s,;])M\s*[:=]\s*(\d+)\b")

# P2 — word forms
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

# P3 — two ints alone on a (mostly) bare comment line; allow leading/trailing
# whitespace and punctuation, but reject if it looks like e.g. atom indices.
_P3 = re.compile(r"^\s*([+-]?\d+)\s+(\d+)\s*$")
# A looser P3 that accepts trailing text like "0 1 RB3LYP/6-31G(d) opt"
_P3_LOOSE = re.compile(r"^\s*([+-]?\d{1,2})\s+([1-9]\d?)\b")

# P4 — S = ...
_P4 = re.compile(r"\bS\s*=\s*(\d+(?:\.\d+)?|\d+/\d+)\b", re.IGNORECASE)


def parse_S_to_mult(s):
    """'0' -> 1, '1/2' -> 2, '1' -> 3, '0.5' -> 2"""
    if "/" in s:
        a, b = s.split("/")
        S = float(a) / float(b)
    else:
        S = float(s)
    mult = int(round(2*S + 1))
    return mult if 1 <= mult <= 10 else None


def parse_comment(line):
    """Try patterns in order. Return (charge, mult, pattern_label) or None."""
    if not line or not line.strip():
        return None
    L = line.strip()
    if len(L) > 500:
        L = L[:500]
    low = L.lower()

    # P1 — explicit fields (charge AND mult independently extractable)
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
    # P1 partial — only one of the two named
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

    # P2 — word forms
    p2_ch = None; p2_mu = None
    # mult word: search whole-word
    for w, val in _MULT_WORDS.items():
        if re.search(r"\b" + w + r"\b", low):
            p2_mu = val
            break
    for w, val in _CHARGE_WORDS.items():
        if re.search(r"\b" + w + r"\b", low):
            p2_ch = val
            break

    # Combine P1 and P2 partials
    if p1_ch is not None or p2_ch is not None:
        ch_candidate = p1_ch if p1_ch is not None else p2_ch
    else:
        ch_candidate = None
    if p1_mu is not None or p2_mu is not None:
        mu_candidate = p1_mu if p1_mu is not None else p2_mu
    else:
        mu_candidate = None

    if ch_candidate is not None and mu_candidate is not None:
        # was P1 or P2 — label by what was used
        label = ("P1_explicit" if (p1_ch is not None and p1_mu is not None)
                 else "P1+P2"  if (p1_ch is not None or p1_mu is not None)
                 else "P2_words")
        return (ch_candidate, mu_candidate, label)

    # P3 — strict "two ints alone"
    m3 = _P3.match(L)
    if m3:
        ch = int(m3.group(1)); mu = int(m3.group(2))
        if -5 <= ch <= 5 and 1 <= mu <= 10:
            return (ch, mu, "P3_two_ints_strict")

    # P4 — spin S =
    m4 = _P4.search(L)
    if m4:
        mu_p4 = parse_S_to_mult(m4.group(1))
        if mu_p4 is not None:
            ch_use = ch_candidate if ch_candidate is not None else 0
            label = ("P4_spin+P1" if p1_ch is not None
                     else "P4_spin+P2" if p2_ch is not None
                     else "P4_spin_only")
            return (ch_use, mu_p4, label)

    # P3 loose — comment line that STARTS with "0 1 ..." pattern
    m3l = _P3_LOOSE.match(L)
    if m3l:
        ch = int(m3l.group(1)); mu = int(m3l.group(2))
        if -5 <= ch <= 5 and 1 <= mu <= 10:
            # require the line to look "header-like": short, no fancy words
            if len(L) < 80 and not any(w in low for w in ("=", "//", "transition", "trans-state", "atom")):
                return (ch, mu, "P3_two_ints_loose")

    # Partial only (just charge or just mult) — still useful as a record
    if ch_candidate is not None and mu_candidate is None:
        return (ch_candidate, None, "partial_charge_only")
    if mu_candidate is not None and ch_candidate is None:
        return (None, mu_candidate, "partial_mult_only")

    return None


def main():
    n_total = 0
    n_blank = 0
    n_matched = 0
    n_partial = 0
    pattern_counter = Counter()
    charge_counter  = Counter()
    mult_counter    = Counter()
    leaf_paper_match = defaultdict(lambda: {"matched": 0, "total": 0})  # by doi
    paper_match_count = defaultdict(lambda: {"matched": 0, "total": 0, "leaves_with_label": set(), "all_leaves": set()})

    # rolling sample of un-matched non-blank comments
    unmatched_sample = []
    UNMATCHED_LIMIT = 500

    # Output CSV streamed
    gt_path = RESULTS / "charge_spin_groundtruth.csv"
    gt_fp = open(gt_path, "w", newline="")
    gt_w = csv.writer(gt_fp)
    gt_w.writerow(["file", "doi", "leaf_dir", "charge", "multiplicity",
                   "pattern", "raw_comment"])

    print("Walking xyz files ...", file=sys.stderr)
    last_print = 0
    for dp, dnames, fnames in os.walk(DOI_ROOT):
        if dnames:  # not a leaf
            continue
        leaf_path = Path(dp)
        doi = doi_from_path(leaf_path)
        if doi is None:
            doi = "?"

        for fname in fnames:
            if not is_target_xyz(fname):
                continue
            fp = leaf_path / fname
            n_total += 1
            paper_match_count[doi]["total"] += 1
            paper_match_count[doi]["all_leaves"].add(str(leaf_path))

            # progress
            if n_total - last_print >= 50000:
                print(f"  ... scanned {n_total:,} files, matched {n_matched:,}, partial {n_partial:,}",
                      file=sys.stderr)
                last_print = n_total

            try:
                with open(fp, errors="replace") as f:
                    f.readline()       # line 1: atom count
                    comment = f.readline()
            except OSError:
                continue

            if comment is None:
                n_blank += 1
                continue
            comment = comment.rstrip("\n\r")
            if not comment.strip():
                n_blank += 1
                continue

            parsed = parse_comment(comment)
            if parsed is None:
                if len(unmatched_sample) < UNMATCHED_LIMIT:
                    unmatched_sample.append((str(fp), comment.strip()[:200]))
                continue

            ch, mu, label = parsed
            pattern_counter[label] += 1
            if label.startswith("partial"):
                n_partial += 1
            else:
                n_matched += 1
                paper_match_count[doi]["matched"] += 1
                paper_match_count[doi]["leaves_with_label"].add(str(leaf_path))
            if ch is not None:
                charge_counter[ch] += 1
            if mu is not None:
                mult_counter[mu] += 1
            gt_w.writerow([str(fp), doi, str(leaf_path),
                           "" if ch is None else ch,
                           "" if mu is None else mu,
                           label,
                           comment.strip()[:200]])
    gt_fp.close()
    print(f"... scanned {n_total:,} files total. Matched {n_matched:,}, partial {n_partial:,}.",
          file=sys.stderr)

    # ── coverage CSV
    cov_path = RESULTS / "charge_spin_coverage.csv"
    with open(cov_path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["doi", "n_xyz_total", "n_xyz_matched",
                    "n_leaves_total", "n_leaves_with_any_label"])
        for doi in sorted(paper_match_count):
            stats = paper_match_count[doi]
            w.writerow([doi, stats["total"], stats["matched"],
                        len(stats["all_leaves"]),
                        len(stats["leaves_with_label"])])

    # ── value distribution CSV
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
        for fp_str, line in unmatched_sample:
            fp.write(f"{fp_str}\n  {line!r}\n\n")

    # ── summary
    n_papers_total       = len(paper_match_count)
    n_papers_with_label  = sum(1 for d in paper_match_count.values() if d["matched"] > 0)
    n_leaves_total       = sum(len(d["all_leaves"]) for d in paper_match_count.values())
    n_leaves_with_label  = sum(len(d["leaves_with_label"]) for d in paper_match_count.values())

    sum_path = RESULTS / "charge_spin_summary.txt"
    with open(sum_path, "w") as fp:
        fp.write("Charge / multiplicity extraction summary\n")
        fp.write("=" * 50 + "\n\n")
        fp.write(f"Total xyz files scanned:                 {n_total:,}\n")
        fp.write(f"  blank or unreadable comment lines:     {n_blank:,}\n")
        fp.write(f"  matched (full charge+mult or P4):      {n_matched:,}  "
                 f"({100*n_matched/n_total:.2f}%)\n")
        fp.write(f"  partial (charge only OR mult only):    {n_partial:,}  "
                 f"({100*n_partial/n_total:.2f}%)\n")
        fp.write(f"  unmatched (non-blank but no signal):   {n_total - n_matched - n_partial - n_blank:,}\n\n")
        fp.write(f"Pattern frequency (full matches only):\n")
        for pat, c in pattern_counter.most_common():
            fp.write(f"  {pat:30s} {c:>9,}\n")
        fp.write(f"\n--- per-paper coverage ---\n")
        fp.write(f"Total papers seen:                   {n_papers_total:,}\n")
        fp.write(f"Papers with at least 1 labeled xyz:  {n_papers_with_label:,}  "
                 f"({100*n_papers_with_label/n_papers_total:.2f}%)\n")
        fp.write(f"Total leaf dirs:                     {n_leaves_total:,}\n")
        fp.write(f"Leaves with at least 1 labeled xyz:  {n_leaves_with_label:,}  "
                 f"({100*n_leaves_with_label/n_leaves_total:.2f}%)\n\n")

        fp.write(f"--- charge distribution (from full + partial matches) ---\n")
        for ch in sorted(charge_counter):
            fp.write(f"  charge={ch:+3d}    {charge_counter[ch]:>9,}\n")
        fp.write(f"\n--- multiplicity distribution ---\n")
        for mu in sorted(mult_counter):
            fp.write(f"  M={mu}    {mult_counter[mu]:>9,}\n")
        fp.write(f"\nOutputs:\n")
        fp.write(f"  {gt_path}\n  {cov_path}\n  {val_path}\n  {um_path}\n  {sum_path}\n")


if __name__ == "__main__":
    main()
