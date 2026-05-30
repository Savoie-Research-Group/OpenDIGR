#!/usr/bin/env python3
import os
import json
import csv
import glob

DOWNLOADS = "/home/z/Desktop/NEW-Wiley/downloads"
OUT_CSV = "/home/z/Desktop/NEW-Wiley/comp_details_summary.csv"


def normalize(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(str(v) for v in value)
    return str(value)


rows = []
subfolders = sorted(os.listdir(DOWNLOADS))

for sub in subfolders:
    cd_dir = os.path.join(DOWNLOADS, sub, "comp_details")
    if not os.path.isdir(cd_dir):
        continue
    json_files = glob.glob(os.path.join(cd_dir, "*-comp_detail.json"))
    if not json_files:
        json_files = glob.glob(os.path.join(cd_dir, "*.json"))
    if not json_files:
        rows.append([sub, "", "", "", "NO_JSON"])
        continue
    jf = json_files[0]
    try:
        with open(jf, "r", encoding="utf-8") as f:
            data = json.load(f)
        entries = data if isinstance(data, list) else [data]
        def collect(field):
            seen, out = set(), []
            for e in entries:
                if not isinstance(e, dict):
                    continue
                v = e.get(field)
                if v is None:
                    continue
                items = v if isinstance(v, list) else [v]
                for it in items:
                    s = str(it).strip()
                    if s and s not in seen:
                        seen.add(s)
                        out.append(s)
            return "; ".join(out)
        lot = collect("level_of_theory")
        basis = collect("basis_set")
        func = collect("functional")
        rows.append([sub, lot, basis, func, ""])
    except Exception as e:
        rows.append([sub, "", "", "", f"ERROR: {e}"])

with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["subfolder", "level_of_theory", "basis_set", "functional", "note"])
    w.writerows(rows)

print(f"Wrote {len(rows)} rows to {OUT_CSV}")
