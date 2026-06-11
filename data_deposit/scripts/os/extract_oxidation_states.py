#!/usr/bin/env python3
"""
Extract transition metal oxidation states from tar.zst reaction archives.

For each archive, reads:
  - IRC_Analysis/finished_first/frame_00000.json  (reactant)
  - IRC_Analysis/finished_last/frame_00000.json   (product)

Oxidation state = el_valence[element] - e
where 'e' is the lone-pair electron count from the JSON.

Transition metals = full d-block (atomic numbers 21-30, 39-48, 57, 72-80).
"""

import csv
import io
import json
import os
import subprocess
import sys
import tarfile
from multiprocessing import Pool

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Valence electrons (from properties.py)
el_valence = {
    'h':1, 'he':2,
    'li':1, 'be':2, 'b':3, 'c':4, 'n':5, 'o':6, 'f':7, 'ne':8,
    'na':1, 'mg':2, 'al':3, 'si':4, 'p':5, 's':6, 'cl':7, 'ar':8,
    'k':1, 'ca':2, 'sc':3, 'ti':4, 'v':5, 'cr':6, 'mn':7, 'fe':8, 'co':9, 'ni':10, 'cu':11, 'zn':12,
    'ga':3, 'ge':4, 'as':5, 'se':6, 'br':7, 'kr':8,
    'rb':1, 'sr':2, 'y':3, 'zr':4, 'nb':5, 'mo':6, 'tc':7, 'ru':8, 'rh':9, 'pd':10, 'ag':11, 'cd':12,
    'in':3, 'sn':4, 'sb':5, 'te':6, 'i':7, 'xe':8,
    'cs':1, 'ba':2, 'la':3, 'hf':4, 'ta':5, 'w':6, 're':7, 'os':8, 'ir':9, 'pt':10, 'au':11, 'hg':12,
    'tl':3, 'pb':4, 'bi':5, 'po':6, 'at':7, 'rn':8,
}

# Full d-block: atomic numbers 21-30, 39-48, 57, 72-80
TRANSITION_METALS = {
    'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn',  # 21-30
    'y', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd',  # 39-48
    'la',                                                           # 57
    'hf', 'ta', 'w', 're', 'os', 'ir', 'pt', 'au', 'hg',         # 72-80
}

REACTANT_PATH = "IRC_Analysis/finished_first/frame_00000.json"
PRODUCT_PATH = "IRC_Analysis/finished_last/frame_00000.json"

BASE_DIR = "/scratch/negishi/li1724/SI-Downloads/SI_Agent/doi_tar_zsts"
OUTPUT_CSV = os.path.join(BASE_DIR, "transition_metal_oxidation_states.csv")

NUM_WORKERS = 64  # plenty for I/O-bound work on 256 cores


def extract_two_members(tar_zst_path):
    """Decompress once, extract both reactant and product JSON bytes."""
    try:
        proc = subprocess.run(
            ["zstd", "-dc", tar_zst_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30,
        )
        if proc.returncode != 0:
            return None, None
        reactant_bytes = None
        product_bytes = None
        with tarfile.open(fileobj=io.BytesIO(proc.stdout), mode="r:") as tf:
            for name, target in [(REACTANT_PATH, "reactant"), (PRODUCT_PATH, "product")]:
                try:
                    f = tf.extractfile(tf.getmember(name))
                    if f is not None:
                        if target == "reactant":
                            reactant_bytes = f.read()
                        else:
                            product_bytes = f.read()
                except KeyError:
                    pass
        return reactant_bytes, product_bytes
    except Exception:
        return None, None


def compute_tm_oxidation_states(atoms_list):
    """Given atoms list from JSON, return string like 'Pt9:4;Fe3:2' for TM atoms."""
    results = []
    for atom in atoms_list:
        el = atom["el"]
        if el.lower() not in TRANSITION_METALS:
            continue
        valence = el_valence.get(el.lower())
        if valence is None:
            continue
        e = atom["e"]
        os_val = valence - e
        atom_id = atom["id"]
        results.append("{}:{}".format(atom_id, os_val))
    return ";".join(results)


def process_one(tar_zst_path):
    """Process a single tar.zst file. Returns a tuple or None."""
    reactant_bytes, product_bytes = extract_two_members(tar_zst_path)

    if reactant_bytes is None and product_bytes is None:
        return None

    reactant_os = "N/A"
    product_os = "N/A"

    if reactant_bytes is not None:
        try:
            data = json.loads(reactant_bytes)
            reactant_os = compute_tm_oxidation_states(data["atoms"])
            if not reactant_os:
                reactant_os = ""
        except Exception:
            reactant_os = "ERROR"

    if product_bytes is not None:
        try:
            data = json.loads(product_bytes)
            product_os = compute_tm_oxidation_states(data["atoms"])
            if not product_os:
                product_os = ""
        except Exception:
            product_os = "ERROR"

    # Skip if no transition metals found on either side (both empty strings)
    if reactant_os == "" and product_os == "":
        return None
    if reactant_os == "N/A" and product_os == "":
        return None
    if reactant_os == "" and product_os == "N/A":
        return None

    return (tar_zst_path, reactant_os, product_os)


def collect_tar_zst_paths():
    """Find all tar.zst files under 10.* DOI directories."""
    paths = []
    for entry in os.listdir(BASE_DIR):
        if not entry.startswith("10."):
            continue
        doi_dir = os.path.join(BASE_DIR, entry)
        if not os.path.isdir(doi_dir):
            continue
        for root, dirs, files in os.walk(doi_dir):
            for fname in files:
                if fname.endswith(".tar.zst"):
                    paths.append(os.path.join(root, fname))
    return paths


def main():
    print("Collecting tar.zst paths...", flush=True)
    paths = collect_tar_zst_paths()
    total = len(paths)
    print(f"Found {total} tar.zst files.", flush=True)

    written = 0
    errors = 0
    skipped = 0

    with open(OUTPUT_CSV, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["tar_zst_path", "reactant_metal_oxidation_states", "product_metal_oxidation_states"])

        with Pool(NUM_WORKERS) as pool:
            for i, result in enumerate(pool.imap_unordered(process_one, paths, chunksize=256), 1):
                if result is None:
                    skipped += 1
                else:
                    writer.writerow(result)
                    written += 1
                if i % 10000 == 0:
                    print(f"  Processed {i}/{total}  (written={written}, skipped={skipped})", flush=True)

    print(f"\nDone. Wrote {written} rows to {OUTPUT_CSV}")
    print(f"Skipped {skipped} (no TMs or missing both JSONs)")


if __name__ == "__main__":
    main()
