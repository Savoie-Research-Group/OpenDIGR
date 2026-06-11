# Charge / spin ground-truth bundle — 2026-05-29

Author-labeled charge/multiplicity ground-truth scraped from xyz comment
lines across the 12-publisher SI corpus (2,275,388 xyz files scanned),
plus the parity-corrected lowest-spin benchmark that revised the headline
coverage from a buggy 71.2% to the correct **86.79%**.

## Layout

```
data/
  charge_spin_groundtruth.csv      (64 MB, ~192K rows: full + partial matches)
  charge_spin_coverage.csv          per-paper coverage roll-up
  charge_spin_value_distribution.csv  charge × multiplicity histogram
  charge_spin_unmatched_sample.txt   audit sample of non-matching comment lines
  parity_corrected_coverage.csv    (32 MB, ~95K rows: per-row electron count + parity)

scripts/
  extract_charge_spin.py           serial version (kept for reference)
  extract_charge_spin_parallel.py  16-worker multiprocessing.Pool / imap_unordered
                                   chunksize=400. This is the version that produced
                                   the data.
  submit_charge_spin.sge           SGE submission wrapper (16 cores, ~34 min)
  parity_corrected_coverage.py     post-processing step. Parses xyz coords, sums Z,
                                   computes N_electrons = sum(Z) - labeled_charge,
                                   M_min = 1 if N_e even else 2. Compares to
                                   author-labeled M. Fixes the 71.2% benchmark bug.

reports/
  charge_spin_summary.txt          scan-level totals + pattern firing breakdown
  parity_corrected_summary.txt     coverage at M_min, +2, +4, ..., per-tier mismatch
  charge_spin.sge.log              SGE wallclock log
  CHARGE_SPIN_GROUNDTRUTH_REPORT.md  full narrative
```

## Methodology in 30 seconds

`extract_charge_spin_parallel.py` scans every xyz file's line 2 (the comment
line) with a tiered regex stack:

  P1 — named fields (charge=N, mult=N, Q=, M=)
  P2 — word forms (neutral/anion/cation, singlet/doublet/.../octet)
  P3 — Gaussian-header "two ints" (strict: whole line; loose: starts-with + word veto)
  P4 — S = x form (S=0, S=1/2 → mult = 2S+1)

Bounds: |charge| ≤ 5, multiplicity ≤ 10. Pool size 16, chunksize 400.

`parity_corrected_coverage.py` then benchmarks the heuristic of "sample
{-1, 0, +1} charges at the parity-correct lowest-spin M_min" against the
ground truth, with the parity correction described above.

## Headline numbers (also in reports/)

- 2,275,388 xyz scanned · 8,872 blank · 96,303 full matches (4.23%) ·
  95,765 partial · 2,958 papers (5.88%) labeled at least one xyz
- Charge sampling {-1, 0, +1} covers **97.3%** of labeled charges
- Lowest-spin **M_min covers 86.79%** (parity-corrected; was incorrectly
  71.2% before the fix). +2 / +4 / +6 / +8 cumulative: 92.5 / 94.4 / 95.0 / 95.3 %
- 4.74% parity mismatch — concentrated in low-confidence tiers
  (P4_spin_only 35%, P3_loose 46%) and nearly absent in strong tiers
  (P3_strict 2.3%, P1 4.3%)

See `reports/CHARGE_SPIN_GROUNDTRUTH_REPORT.md` for full discussion.
