#!/usr/bin/env python3
"""
Extract DFT single-point energies and xTB IRC barriers from doi_zips/ archives.

Usage:
    python extract_energies.py --keys keys.txt
                               --zip-dir doi_zips/
                               -o energies.csv
                               [-j 16]

For each DOI key (e.g. "10.1039/C1CC15737J/c1cc15737j/06_-1_2"):
  - Opens <zip-dir>/<key>.zip
  - Reads 4 ORCA .out files from <reaction>/DFT-SinglePoint/:
        input.out, finished_first.out, finished_last.out, ts_final_geometry.out
    Extracts "FINAL SINGLE POINT ENERGY  <value>" (Hartree)
  - Reads any <reaction>/*-irc.out (xTB IRC) and extracts Left/TS/Right
    "<X>:  <value> kJ mol" (kJ/mol)

Writes one CSV row per key. Missing values are blank.
"""

import os
import re
import sys
import csv
import zipfile
import argparse
from multiprocessing import Pool

DFT_CLASSES = (
    "input",
    "finished_first",
    "finished_first_opt",
    "finished_last",
    "finished_last_opt",
    "ts_final_geometry",
)

KJ_PER_KCAL = 4.184
HARTREE_TO_KCAL = 627.5094740631

RE_FINAL_SP = re.compile(r"^\s*FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)\s*$")
RE_IRC = re.compile(
    r"^\s*(Left|TS|Right)\s*:\s*(-?\d+(?:\.\d+)?)\s+kJ\s*mol", re.IGNORECASE
)


def detect_reaction_dir(namelist):
    """Return the single top-level directory inside the zip (e.g. '06_-1_2/')."""
    tops = {n.split("/", 1)[0] for n in namelist if n}
    if len(tops) == 1:
        return next(iter(tops)) + "/"
    # fallback: pick the most common
    from collections import Counter
    c = Counter(n.split("/", 1)[0] for n in namelist if "/" in n)
    return (c.most_common(1)[0][0] + "/") if c else ""


def parse_orca_final_sp(text):
    """Return the FINAL SINGLE POINT ENERGY (last occurrence) in Hartree, or None."""
    last = None
    for line in text.splitlines():
        m = RE_FINAL_SP.match(line)
        if m:
            last = float(m.group(1))
    return last


def parse_irc_block(text):
    """Return dict {'Left': float, 'TS': float, 'Right': float} in kJ/mol, or {}."""
    found = {}
    for line in text.splitlines():
        m = RE_IRC.match(line)
        if m:
            key = m.group(1).capitalize() if m.group(1).upper() != "TS" else "TS"
            # Take the *last* hit per key (final summary block at end of file)
            found[key] = float(m.group(2))
    return found


def extract_one(args):
    key, zip_path = args
    row = {"key": key}
    for c in DFT_CLASSES:
        row[f"E_{c}_Eh"] = ""
    row["irc_file"] = ""
    row["IRC_Left_kJmol"] = ""
    row["IRC_TS_kJmol"] = ""
    row["IRC_Right_kJmol"] = ""
    row["IRC_Left_kcalmol"] = ""
    row["IRC_TS_kcalmol"] = ""
    row["IRC_Right_kcalmol"] = ""
    row["DFT_barrier_from_first_kcalmol"] = ""
    row["DFT_barrier_from_first_opt_kcalmol"] = ""
    row["DFT_barrier_from_last_kcalmol"] = ""
    row["DFT_barrier_from_last_opt_kcalmol"] = ""
    row["error"] = ""

    try:
        if not os.path.exists(zip_path):
            row["error"] = "zip_missing"
            return row

        with zipfile.ZipFile(zip_path, "r") as zf:
            names = zf.namelist()
            rxn = detect_reaction_dir(names)

            # DFT single points
            for c in DFT_CLASSES:
                target = f"{rxn}DFT-SinglePoint/{c}.out"
                if target in zf.namelist():
                    with zf.open(target) as f:
                        text = f.read().decode("utf-8", errors="replace")
                    e = parse_orca_final_sp(text)
                    if e is not None:
                        row[f"E_{c}_Eh"] = f"{e:.10f}"

            # xTB IRC: any top-level <reaction>/*-irc.out (skip DFT-SinglePoint/IRC_Analysis)
            irc_candidates = [
                n for n in names
                if n.startswith(rxn) and n.endswith("-irc.out")
                and "/DFT-SinglePoint/" not in n
                and "/IRC_Analysis/" not in n
                and n.count("/") == rxn.count("/")  # at the top level inside rxn/
            ]
            if irc_candidates:
                irc_name = irc_candidates[0]
                row["irc_file"] = os.path.basename(irc_name)
                with zf.open(irc_name) as f:
                    text = f.read().decode("utf-8", errors="replace")
                d = parse_irc_block(text)
                row["IRC_Left_kJmol"]  = f"{d['Left']:.2f}"  if "Left"  in d else ""
                row["IRC_TS_kJmol"]    = f"{d['TS']:.2f}"    if "TS"    in d else ""
                row["IRC_Right_kJmol"] = f"{d['Right']:.2f}" if "Right" in d else ""
                if "Left"  in d: row["IRC_Left_kcalmol"]  = f"{d['Left']  / KJ_PER_KCAL:.4f}"
                if "TS"    in d: row["IRC_TS_kcalmol"]    = f"{d['TS']    / KJ_PER_KCAL:.4f}"
                if "Right" in d: row["IRC_Right_kcalmol"] = f"{d['Right'] / KJ_PER_KCAL:.4f}"

        # DFT barriers in kcal/mol: TS - endpoint, using both the IRC-geometry SP
        # (finished_first/last) and the DFT-reoptimized endpoint (finished_first_opt/last_opt).
        # The *_opt versions are the apples-to-apples comparison with the xTB IRC barriers.
        e_ts        = row.get("E_ts_final_geometry_Eh") or ""
        e_first     = row.get("E_finished_first_Eh") or ""
        e_first_opt = row.get("E_finished_first_opt_Eh") or ""
        e_last      = row.get("E_finished_last_Eh") or ""
        e_last_opt  = row.get("E_finished_last_opt_Eh") or ""
        if e_ts and e_first:
            row["DFT_barrier_from_first_kcalmol"]     = f"{(float(e_ts) - float(e_first))     * HARTREE_TO_KCAL:.4f}"
        if e_ts and e_first_opt:
            row["DFT_barrier_from_first_opt_kcalmol"] = f"{(float(e_ts) - float(e_first_opt)) * HARTREE_TO_KCAL:.4f}"
        if e_ts and e_last:
            row["DFT_barrier_from_last_kcalmol"]      = f"{(float(e_ts) - float(e_last))      * HARTREE_TO_KCAL:.4f}"
        if e_ts and e_last_opt:
            row["DFT_barrier_from_last_opt_kcalmol"]  = f"{(float(e_ts) - float(e_last_opt))  * HARTREE_TO_KCAL:.4f}"

        return row
    except Exception as e:
        row["error"] = type(e).__name__ + ": " + str(e)
        return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keys", required=True, help="File with one DOI key per line")
    ap.add_argument("--zip-dir", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()

    with open(args.keys) as f:
        keys = [ln.strip() for ln in f if ln.strip()]
    pairs = [(k, os.path.join(args.zip_dir, k + ".zip")) for k in keys]
    print(f"Extracting {len(pairs)} zips with {args.j} workers ...")

    fieldnames = [
        "key",
        "E_input_Eh",
        "E_finished_first_Eh",
        "E_finished_first_opt_Eh",
        "E_finished_last_Eh",
        "E_finished_last_opt_Eh",
        "E_ts_final_geometry_Eh",
        "DFT_barrier_from_first_kcalmol",
        "DFT_barrier_from_first_opt_kcalmol",
        "DFT_barrier_from_last_kcalmol",
        "DFT_barrier_from_last_opt_kcalmol",
        "irc_file",
        "IRC_Left_kJmol",
        "IRC_TS_kJmol",
        "IRC_Right_kJmol",
        "IRC_Left_kcalmol",
        "IRC_TS_kcalmol",
        "IRC_Right_kcalmol",
        "error",
    ]

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames)
        w.writeheader()

        if args.j == 1:
            for i, p in enumerate(pairs, 1):
                w.writerow(extract_one(p))
                if i % 1000 == 0:
                    print(f"  {i}/{len(pairs)}")
        else:
            with Pool(args.j) as pool:
                for i, row in enumerate(pool.imap_unordered(extract_one, pairs, chunksize=32), 1):
                    w.writerow(row)
                    if i % 1000 == 0:
                        print(f"  {i}/{len(pairs)}")

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
