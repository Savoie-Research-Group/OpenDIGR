# Charge / multiplicity ground-truth extraction from xyz comment lines

## Why this exists

A reviewer questioned the reliability of our automated charge/spin
assignment pipeline (the "two-pass" approach that samples charges
{−1, 0, +1} with lowest-spin multiplicity, then explores higher
spins via a WBO scan). To respond, we built a ground-truth set from
**author-labeled** charge/multiplicity annotations in xyz comment
lines across the entire DIGR corpus.

This report describes the extraction methodology and what we found.
The companion CSV (`charge_spin_groundtruth.csv`) is the artifact you
can cross-join against the automated assignments to compute a true
accuracy / confusion matrix.

## What was scanned

| Metric | Count |
|---|---:|
| `.xyz` files scanned (excluding `*repacked.xyz`) | 2,275,388 |
| Leaf directories walked | 261,613 |
| DOIs with ≥ 1 xyz file | 50,338 |
| Wallclock on 16-core SGE node | 34 minutes |

The scan covered every xyz file under
`/groups/bsavoie2/zli43/Original-Files/doi_files/<DOI>/{pdf-processed-XYZ,Raw-XYZ}/.../leaf_dir/`
(both publishers' real coordinates and PDF-extracted ones).

## Extraction strategy

For each xyz file we read **only line 2** — the standard xyz
"comment" line that authors may use to annotate the structure.

Standard xyz format:
```
<n_atoms>          ← line 1, ignored here
<comment>          ← line 2, parsed
<element> <x> <y> <z>
...
```

### Tiered regex stack

Patterns are applied in descending order of confidence. The first
tier that produces both a charge and a multiplicity wins; if a tier
matches only one of the two, we fall through to lower tiers and
combine partials when possible.

| tier | what it matches | example comment |
|---|---|---|
| **P1 — explicit named fields** | `charge[:=]?<int>` together with `mult[:=]?<int>` (also `q=`, `chg=`, `multiplicity=`, `2s+1=`, `M=<int>`) | `charge=0 mult=1, TS for Diels-Alder` |
| **P2 — word forms** | charge words (`neutral`/`anion`/`cation`/`dianion`/...) → 0/−1/+1/−2; mult words (`singlet`/`doublet`/...`/octet`) → 1/2/.../8 | `neutral singlet, optimized at B3LYP/6-31G(d)` |
| **P1 + P2 mix** | combine partials from the two | `charge=−1, triplet` |
| **P3 — two ints strict** | whole comment is exactly `[+-]?<int> <int>` (Gaussian-header convention) | `0 1` |
| **P4 — spin S = …** | `S = 0`, `S = 1/2`, `S = 1` → multiplicity = 2S+1 (charge taken from P1/P2 partial if present, else assumed 0) | `S=1/2` |
| **P3 loose** | comment *starts* with two small ints and looks header-like (short, no `=`, no chemistry words) | `0 1 RB3LYP/6-31G(d) opt` |
| **partial only** | only charge OR only multiplicity extractable | `cationic`, `multiplicity 3` |

Sanity bounds enforced after every match: charge ∈ [−5, +5],
multiplicity ∈ [1, 10]. Anything outside those bounds is rejected as a
false-positive (e.g. atom indices that happened to fit `<int> <int>`).

## Coverage

### File / paper / leaf-dir level

| level | total | labeled | rate |
|---|---:|---:|---:|
| files (full charge+mult match) | 2,275,388 | **96,303** | **4.23%** |
| files (charge or mult, partial OK) | 2,275,388 | 192,068 | 8.44% |
| leaf directories (≥ 1 labeled file) | 261,613 | **38,331** | **14.65%** |
| papers / DOIs (≥ 1 labeled file) | 50,338 | **2,958** | **5.88%** |
| comment lines blank or unreadable | — | 8,872 | 0.39% |

**Headline:** roughly **6 % of papers** in the DIGR corpus
explicitly annotate charge or multiplicity in their xyz comments.
The bulk are silent.

### Pattern firings (full matches only)

| pattern | count | % of matches |
|---|---:|---:|
| **P3_two_ints_strict** (`0 1`) | 73,333 | 76.2% |
| **P1_explicit** (`charge=0 mult=1`) | 14,583 | 15.1% |
| P4_spin_only (`S=…`) | 4,502 | 4.7% |
| P3_two_ints_loose (`0 1 RB3LYP/...`) | 1,746 | 1.8% |
| P4_spin+P1 + P4_spin+P2 | 957 | 1.0% |
| P2_words (`neutral singlet`) | 733 | 0.8% |
| P1+P2 (mixed) | 449 | 0.5% |

**Major finding:** the dominant explicit convention is **not**
`charge=… mult=…` — it's the bare **`<int> <int>`** Gaussian-header
form. Authors paste the Gaussian input header directly into the xyz
comment far more often than they write out the named fields.

### Partial-only matches (one of two)

| pattern | count |
|---|---:|
| partial_charge_only (e.g. just `cationic`) | 78,147 |
| partial_mult_only (e.g. just `multiplicity 3`) | 17,618 |

These partials are useful as labeled charge data or labeled spin
data, but they don't constitute full ground-truth rows in
`charge_spin_groundtruth.csv` (their `charge` or `multiplicity`
column is blank). For the strict benchmark, use only the 96,303 full
matches.

### How concentrated is the labeled data?

Among the 2,958 labeled papers:

| files-labeled per paper | papers | fraction |
|---|---:|---:|
| 1 file | 558 | 18.9% |
| 2–5 files | 715 | 24.2% |
| 6–25 files | 1,005 | 34.0% |
| 26–100 files | 556 | 18.8% |
| > 100 files | 124 | 4.2% |

Labeled papers tend to label many of their structures, not just one
— a quarter have 6+ labeled xyz, and 124 papers contribute 100+
each. Treating labels-per-paper as independent observations
overstates the breadth: ~4% of papers contribute disproportionately
to the file-level count.

For unbiased per-paper benchmarking, dedupe to one row per
`(DOI, charge, multiplicity)` triple first.

## Value distributions (from labeled set)

### Charge

| charge | count | % of charge-labeled files |
|---:|---:|---:|
| −5 | 7 | 0.0% |
| −4 | 75 | 0.0% |
| −3 | 264 | 0.2% |
| −2 | 883 | 0.5% |
| **−1** | 8,022 | 4.6% |
| **0** | 96,257 | 55.1% |
| **+1** | 65,552 | 37.5% |
| +2 | 2,426 | 1.4% |
| +3 | 469 | 0.3% |
| +4 | 288 | 0.2% |
| +5 | 207 | 0.1% |
| **{−1, 0, +1} combined** | **169,831** | **97.3%** |

### Multiplicity

| multiplicity | count | % of mult-labeled files |
|---:|---:|---:|
| **M = 1 (singlet)** | 81,052 | **71.2%** |
| M = 2 (doublet) | 13,221 | 11.6% |
| M = 3 (triplet) | 11,355 | 10.0% |
| M = 4 (quartet) | 2,279 | 2.0% |
| M = 5 (quintet) | 3,311 | 2.9% |
| M = 6 (sextet) | 1,000 | 0.9% |
| M = 7 (septet) | 1,100 | 1.0% |
| M = 8 | 217 | 0.2% |
| M = 9 | 267 | 0.2% |
| M = 10 | 119 | 0.1% |

### Top 10 (charge, multiplicity) combinations (full matches only)

| (q, M) | description | count |
|---|---|---:|
| (0, 1) | neutral singlet | 60,938 |
| (0, 2) | neutral doublet | 7,428 |
| (+1, 1) | cation singlet | 7,256 |
| (−1, 1) | anion singlet | 5,228 |
| (0, 3) | neutral triplet | 3,911 |
| (+1, 2) | cation doublet | 1,757 |
| (0, 5) | neutral quintet | 1,260 |
| (−1, 2) | anion doublet | 841 |
| (0, 7) | neutral septet | 706 |
| (0, 4) | neutral quartet | 659 |
| (+2, 1) | dication singlet | 620 |
| (+1, 3) | cation triplet | 606 |

Closed-shell neutrals dominate, as expected for organic
mechanisms. The high-spin tail (M = 5/7 with q = 0) is small but
not negligible — likely first-row transition-metal complexes.

## Implications for the automated assignment

Two specific implications for your `{−1, 0, +1} × lowest-spin`
sampling, then `WBO` higher-spin pass:

### Charge sampling — well-calibrated

**97.3 % of all author-labeled charges fall in {−1, 0, +1}.** The
remaining 2.7 % are mostly ±2 with a small ±3 / ±4 / ±5 tail.
Coverage by sampling depth:

| sampling depth | coverage |
|---|---:|
| {0} only (neutrals) | 55 % |
| {−1, 0, +1} | **97.3 %** |
| {−2, −1, 0, +1, +2} | **99.2 %** |
| {−3 … +3} | 99.7 % |

Your current 3-charge sampling captures all but 2.7 % of cases the
authors themselves bothered to label. Expanding to ±2 would catch
another ~1.9 percentage points (likely the polynuclear
transition-metal cluster papers in coordination chemistry).

### Multiplicity sampling — "lowest-spin first" misses 28.8 %

| sampling depth | coverage of labeled multiplicities |
|---|---:|
| M = 1 only | 71.2 % |
| M ∈ {1, 2} (singlet + doublet) | 82.8 % |
| M ∈ {1, 2, 3} | **92.7 %** |
| M ∈ {1, …, 5} | 97.6 % |
| M ∈ {1, …, 7} | 99.5 % |

The two-pass WBO scan is exactly the mitigation reviewers will ask
about for the **28.8 %** of cases where `M = 1` is wrong. The
benchmark set lets you quantify how often the WBO scan recovers the
correct spin state for those cases — that's the headline number for
the reviewer response.

The M = 5 / M = 7 high-spin populations (4,411 files combined)
deserve a dedicated drill-down by element class — these are almost
certainly Fe(II/III), Mn(II/III), and Cr(III) sites where the
"lowest-spin default" is systematically incorrect.

## How to use `charge_spin_groundtruth.csv`

One row per xyz file that matched any pattern. Columns:

```
file                full path to the xyz
doi                 paper DOI (10.XXXX/yyyy)
leaf_dir            parent directory (= one "molecule" in our taxonomy)
charge              integer in [−5, +5], or blank for partial_mult_only
multiplicity        integer in [1, 10], or blank for partial_charge_only
pattern             which extraction tier fired (P1_explicit, P3_two_ints_strict, …)
raw_comment         first 200 chars of the line-2 comment (provenance)
```

### Recommended benchmark workflow

1. **Filter to full matches**: `pattern NOT LIKE 'partial%'`. That's
   your 96,303-row ground-truth set.
2. **Optionally restrict to highest-confidence patterns**:
   `pattern IN ('P1_explicit', 'P1+P2', 'P2_words', 'P4_spin+P1', 'P4_spin+P2')`.
   That's ~16,500 rows where the author *named* both fields. P3
   matches are still trustworthy in aggregate, but each P3 row has
   slightly higher false-positive risk than a P1 row.
3. **Dedupe**: collapse to one row per `(doi, leaf_dir, charge,
   multiplicity)` so heavy-labeling papers don't dominate.
4. **Join with the automated assignments**: use `(doi, leaf_dir)`
   (or the full file path if available) as the key.
5. **Compute** per-(q, M) accuracy plus a charge × multiplicity
   confusion matrix.

## Caveats

- **The 6 % labeled subset is not a random sample of the full
  corpus.** Authors who annotate xyz comments may systematically
  differ from those who don't (e.g. authors using ORCA / xTB / open
  workflows tend to embed Gaussian-style headers more often). The
  absolute distributions reported here reflect that subset, not
  the corpus as a whole.

- **P3 "two ints" matches are a convention bet.** A bare `0 1` line
  almost always means charge 0, multiplicity 1, but in principle
  it could be unrelated integers (atom indices, restraint IDs,
  ...). Sanity bounds (charge ≤ |5|, mult ≤ 10) reject the most
  egregious cases, but a few false positives likely survive. P1
  ("`charge=… mult=…`") is the safest tier.

- **We cannot independently verify the author's labels.** If an
  author wrote `charge=0` for a system that is actually +1, our
  ground truth inherits their mistake. The right framing of the
  benchmark is "agreement with author-claimed assignments" rather
  than "agreement with truth".

- **Partial matches (charge-only or multiplicity-only) carry less
  information.** They're in the CSV for completeness but should
  not be used for full (q, M) accuracy scoring.

- **The unlabeled 94 % of papers are not characterized here.** The
  reviewer may push back: are silent papers more likely to be
  organic (closed-shell, neutral) or more likely to include exotic
  high-spin chemistry? Investigating that would require a separate
  spot-check exercise — pick N silent papers, look at the actual
  chemistry described, and tally.

## Files

| File | Contents |
|---|---|
| `extract_charge_spin_parallel.py` | The parallel extraction script (16-worker `Pool.imap_unordered`). |
| `extract_charge_spin.py` | Original single-threaded version (kept for reference). |
| `submit_charge_spin.sge` | SGE submission wrapper (`-pe smp 16`, 2h walltime). |
| **`charge_spin_groundtruth.csv`** | **One row per labeled xyz — the benchmark set.** 192K rows, 65 MB. |
| `charge_spin_coverage.csv` | Per-paper: total xyz / matched xyz / total leaves / leaves-with-any-label. |
| `charge_spin_value_distribution.csv` | Long-form: charge histogram, multiplicity histogram, pattern-firing counts. |
| `charge_spin_unmatched_sample.txt` | First 500 non-blank comments that did NOT match — for tuning. |
| `charge_spin_summary.txt` | Auto-generated text summary (counts + percentages). |
| `charge_spin.sge.log` | SGE job log from the run that produced the above. |
| **`CHARGE_SPIN_GROUNDTRUTH_REPORT.md`** | **This document.** |

## Reproduce

```bash
cd /groups/bsavoie2/zli43/GoldDIGR-Comp-details/results

# Re-run on the cluster (16 cores, ~34 min)
qsub submit_charge_spin.sge

# Or locally with N workers (small N if testing)
python3 extract_charge_spin_parallel.py 4
```

To extend the regex stack: edit the `_P1_*`, `_MULT_WORDS`,
`_CHARGE_WORDS`, `_P3*`, `_P4` definitions at the top of
`extract_charge_spin_parallel.py` and re-run. The pattern stack is
applied in code order; higher-confidence tiers should come first.
