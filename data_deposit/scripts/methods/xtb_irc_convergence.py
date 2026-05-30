#!/usr/bin/env python3
"""xTB-IRC convergence rate, binned by original DFT functional and element class.

Definitions
-----------
Reaction (denominator):
    A unique reaction key appearing in `All-reactions.txt`. That file
    lists every TS for which the xTB-IRC produced bond changes ("a
    reaction happened"). 472,120 keys.

Converged (numerator):
    The same reaction key also appears in
    `052026-all-DFT-IRC-energies-cleaned.csv`, where every barrier
    column (IRC_Left/TS/Right + four DFT barriers) is populated. All
    rows in the CSV are populated (no missing values in the cleaned
    file), so "in the CSV" ≡ "fully converged". 84,355 keys overlap.

Per-paper attributes (from paper_methods.csv):
    Each DOI is annotated with the canonicalized list of DFT
    functional(s) the original paper used and the set of elements
    present in its xyz files. A reaction inherits its paper's
    attributes. Papers list multiple functionals get one reaction
    credit per functional (and similarly for element classes).

Outputs (in results/)
---------------------
    xtb_irc_convergence_per_functional.csv
    xtb_irc_convergence_per_basis.csv
    xtb_irc_convergence_per_element_class.csv
    xtb_irc_convergence_per_software.csv
    xtb_irc_convergence_per_functional.png   bar chart
    xtb_irc_convergence_per_element_class.png  bar chart
    xtb_irc_convergence_summary.txt          overall figures
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")
OTHER   = RESULTS / "Other-files"
CSV_PATH    = OTHER / "052026-all-DFT-IRC-energies-cleaned.csv"
RX_PATH     = OTHER / "All-reactions.txt"
METHODS_CSV = RESULTS / "paper_methods.csv"

# Minimum N (reactions in the denominator) to include a bin in the per-bin tables.
MIN_N = 200


def doi_from_key(key):
    parts = key.split("/", 2)
    if len(parts) >= 2:
        return parts[0] + "/" + parts[1]
    return None


def parse_methods():
    """{doi -> dict with set('functionals','basis_sets','softwares','elements')}"""
    out = {}
    with open(METHODS_CSV) as fp:
        r = csv.DictReader(fp)
        for row in r:
            out[row["doi"]] = {
                "functionals": set(s for s in row["functionals"].split(" | ") if s),
                "basis_sets":  set(s for s in row["basis_sets"].split(" | ") if s),
                "softwares":   set(s for s in row["softwares"].split(" | ") if s),
                "elements":    set(row["elements"].split()),
            }
    return out


# Element class
TM_3D = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn".split())
TM_4D = set("Y Zr Nb Mo Tc Ru Rh Pd Ag Cd".split())
TM_5D = set("Hf Ta W Re Os Ir Pt Au Hg".split())
LANTHANIDES = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu".split())
ACTINIDES   = set("Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr".split())

CLASS_ORDER = ["main_group_only", "3d_TM", "4d_TM", "5d_TM",
               "lanthanides", "actinides", "no_xyz_data"]

def element_class(els):
    if not els:                return "no_xyz_data"
    if els & ACTINIDES:        return "actinides"
    if els & LANTHANIDES:      return "lanthanides"
    if els & TM_5D:            return "5d_TM"
    if els & TM_4D:            return "4d_TM"
    if els & TM_3D:            return "3d_TM"
    return "main_group_only"


def write_rate_csv(path, header_label, counts_total, counts_converged):
    rows = []
    for k, tot in counts_total.items():
        conv = counts_converged.get(k, 0)
        rate = conv / tot if tot else 0.0
        rows.append((k, tot, conv, rate))
    rows.sort(key=lambda r: -r[1])
    with open(path, "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow([header_label, "n_reactions", "n_converged", "convergence_rate"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], f"{r[3]:.4f}"])
    return rows


def bar_chart(rows, title, xlabel, outpath, color="#4477AA", top_n=15,
              show_pct_threshold_n=None):
    """rows = (label, total, converged, rate). Bars are convergence rate;
    annotate with n/total."""
    rows = [r for r in rows if r[1] >= (show_pct_threshold_n or 0)]
    rows = rows[:top_n]
    if not rows:
        return
    labels = [f"{r[0]}\n(n={r[1]:,})" for r in rows]
    rates  = [r[3]*100 for r in rows]
    fig, ax = plt.subplots(figsize=(11, max(5, len(rows)*0.35 + 1)))
    bars = ax.barh(range(len(labels)), rates, color=color,
                   edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    xmax = max(rates) if rates else 1
    for i, (rate, n, c) in enumerate(zip(rates, [r[1] for r in rows],
                                          [r[2] for r in rows])):
        ax.text(rate + xmax*0.005, i, f"{rate:.1f}% ({c:,}/{n:,})",
                va="center", ha="left", fontsize=9)
    ax.set_xlim(0, max(xmax*1.18, 30))
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)
    print(f"wrote {outpath}")


def main():
    print("Loading paper_methods.csv ...", file=sys.stderr)
    methods = parse_methods()
    print(f"  {len(methods)} DOIs", file=sys.stderr)

    print("Loading All-reactions.txt ...", file=sys.stderr)
    rx_keys = set()
    with open(RX_PATH) as fp:
        for line in fp:
            if line.startswith("#") or not line.strip():
                continue
            rx_keys.add(line.split("\t")[0])
    print(f"  {len(rx_keys)} denominator keys", file=sys.stderr)

    print("Loading converged keys from CSV ...", file=sys.stderr)
    conv_keys = set()
    with open(CSV_PATH) as fp:
        r = csv.DictReader(fp)
        for row in r:
            conv_keys.add(row["key"])
    print(f"  {len(conv_keys)} converged keys in CSV", file=sys.stderr)

    # Subset of converged keys that are also "reactions" (denominator).
    converged_and_reacted = conv_keys & rx_keys
    print(f"  overlap (converged AND in All-reactions.txt): "
          f"{len(converged_and_reacted)}", file=sys.stderr)

    # ── overall convergence rate
    overall_total = len(rx_keys)
    overall_conv  = len(converged_and_reacted)
    overall_rate  = overall_conv / overall_total if overall_total else 0

    # ── per-functional / basis / software / element-class bins
    fn_total = defaultdict(int); fn_conv = defaultdict(int)
    bs_total = defaultdict(int); bs_conv = defaultdict(int)
    sw_total = defaultdict(int); sw_conv = defaultdict(int)
    cls_total = defaultdict(int); cls_conv = defaultdict(int)
    fnxcls_total = defaultdict(int); fnxcls_conv = defaultdict(int)

    n_no_doi_match = 0

    for key in rx_keys:
        doi = doi_from_key(key)
        is_conv = key in conv_keys
        paper = methods.get(doi)
        if paper is None:
            n_no_doi_match += 1
            # Still credit the element-class bucket as no_xyz_data so we
            # don't drop the row entirely from the class breakdown.
            cls_total["no_paper_metadata"] += 1
            if is_conv:
                cls_conv["no_paper_metadata"] += 1
            continue

        functionals = paper["functionals"] or {"(no functional listed)"}
        basis_sets  = paper["basis_sets"]  or {"(no basis listed)"}
        softwares   = paper["softwares"]   or {"(no software listed)"}
        cls         = element_class(paper["elements"])

        for fn in functionals:
            fn_total[fn] += 1
            if is_conv: fn_conv[fn] += 1
            fnxcls_total[(fn, cls)] += 1
            if is_conv: fnxcls_conv[(fn, cls)] += 1
        for bs in basis_sets:
            bs_total[bs] += 1
            if is_conv: bs_conv[bs] += 1
        for sw in softwares:
            sw_total[sw] += 1
            if is_conv: sw_conv[sw] += 1
        cls_total[cls] += 1
        if is_conv: cls_conv[cls] += 1

    print(f"\nreactions whose DOI had no paper_methods row: {n_no_doi_match}",
          file=sys.stderr)

    # ── write tables
    fn_rows  = write_rate_csv(RESULTS / "xtb_irc_convergence_per_functional.csv",
                              "functional", fn_total, fn_conv)
    bs_rows  = write_rate_csv(RESULTS / "xtb_irc_convergence_per_basis.csv",
                              "basis_set", bs_total, bs_conv)
    sw_rows  = write_rate_csv(RESULTS / "xtb_irc_convergence_per_software.csv",
                              "software", sw_total, sw_conv)
    cls_rows = write_rate_csv(RESULTS / "xtb_irc_convergence_per_element_class.csv",
                              "element_class", cls_total, cls_conv)

    # functional × element class joint
    rows = []
    for (fn, cls), tot in fnxcls_total.items():
        conv = fnxcls_conv.get((fn, cls), 0)
        rows.append((fn, cls, tot, conv, conv/tot if tot else 0.0))
    rows.sort(key=lambda r: -r[2])
    with open(RESULTS / "xtb_irc_convergence_functional_x_class.csv",
              "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["functional", "element_class", "n_reactions",
                    "n_converged", "convergence_rate"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.4f}"])

    # ── plots
    # filter to bins with enough N to be statistically meaningful;
    # also drop the "(no <field> listed)" placeholders.
    def _named(rows):
        return [r for r in rows
                if r[1] >= MIN_N and not r[0].startswith("(no ")]
    fn_plot  = _named(fn_rows)
    bs_plot  = _named(bs_rows)
    sw_plot  = _named(sw_rows)
    cls_plot = [r for r in cls_rows if r[1] >= MIN_N]

    bar_chart(fn_plot,
              f"xTB-IRC convergence rate by original DFT functional\n"
              f"(overall = {overall_rate*100:.1f}%, N_total = {overall_total:,}; "
              f"papers without a listed functional excluded)",
              "convergence rate (%)",
              RESULTS / "xtb_irc_convergence_per_functional.png", "#4477AA")

    bar_chart(bs_plot,
              f"xTB-IRC convergence rate by original DFT basis set\n"
              f"(overall = {overall_rate*100:.1f}%; papers without a listed basis excluded)",
              "convergence rate (%)",
              RESULTS / "xtb_irc_convergence_per_basis.png", "#117733")

    # element-class plot: exclude the two "non-chemistry" buckets
    cls_plot_clean = [r for r in cls_plot
                      if r[0] not in ("no_xyz_data", "no_paper_metadata")]
    bar_chart(cls_plot_clean,
              f"xTB-IRC convergence rate by element class\n"
              f"(overall = {overall_rate*100:.1f}%; "
              f"unclassified buckets excluded)",
              "convergence rate (%)",
              RESULTS / "xtb_irc_convergence_per_element_class.png", "#CC6677",
              top_n=20)

    bar_chart(sw_plot,
              f"xTB-IRC convergence rate by original software\n"
              f"(overall = {overall_rate*100:.1f}%; papers without a listed software excluded)",
              "convergence rate (%)",
              RESULTS / "xtb_irc_convergence_per_software.png", "#882255")

    # ── summary file
    with open(RESULTS / "xtb_irc_convergence_summary.txt", "w") as fp:
        fp.write("xTB-IRC convergence summary\n")
        fp.write("=" * 40 + "\n\n")
        fp.write(f"Denominator = unique reaction keys in All-reactions.txt:  {overall_total:,}\n")
        fp.write(f"Numerator   = same keys appearing in cleaned-DFT-IRC CSV: {overall_conv:,}\n")
        fp.write(f"Overall convergence rate:  {overall_rate*100:.2f}%\n\n")

        fp.write(f"Reactions whose DOI had no paper_methods row: {n_no_doi_match:,}\n")
        fp.write(f"(These were counted under 'no_paper_metadata' "
                 "in the element-class table.)\n\n")

        fp.write("Top 15 functionals by N (with their convergence rates):\n")
        for r in fn_rows[:15]:
            fp.write(f"  {r[0]:<28s} N={r[1]:>7,}  conv={r[2]:>7,}  rate={r[3]*100:5.1f}%\n")
        fp.write("\nElement-class breakdown:\n")
        for r in cls_rows:
            fp.write(f"  {r[0]:<28s} N={r[1]:>7,}  conv={r[2]:>7,}  rate={r[3]*100:5.1f}%\n")
        fp.write("\nNotes:\n")
        fp.write("- A paper that lists multiple functionals contributes its\n"
                 "  reactions to EACH listed functional's bin (so per-functional\n"
                 "  N values can sum to more than 'overall_total').\n")
        fp.write("- 'no_xyz_data' = paper had a comp_details row but its xyz\n"
                 "  files could not be parsed (or had no DOI match).\n")


if __name__ == "__main__":
    main()
