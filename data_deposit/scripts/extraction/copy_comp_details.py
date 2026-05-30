#!/usr/bin/env python3
import os
import shutil

SOURCES = {
    "New-Wiley": "/home/z/Desktop/NEW-Wiley/downloads",
    "Wiley": "/home/z/Desktop/Wiley/downloads",
}
DEST_ROOT = "/home/z/Desktop/Comp_details"

os.makedirs(DEST_ROOT, exist_ok=True)

for label, src_root in SOURCES.items():
    dest_label_root = os.path.join(DEST_ROOT, label)
    os.makedirs(dest_label_root, exist_ok=True)
    copied = 0
    skipped = 0
    subs = os.listdir(src_root)
    for sub in subs:
        src_cd = os.path.join(src_root, sub, "comp_details")
        if not os.path.isdir(src_cd):
            skipped += 1
            continue
        dest_sub = os.path.join(dest_label_root, sub)
        os.makedirs(dest_sub, exist_ok=True)
        dest_cd = os.path.join(dest_sub, "comp_details")
        if os.path.exists(dest_cd):
            shutil.rmtree(dest_cd)
        shutil.copytree(src_cd, dest_cd)
        copied += 1
    print(f"{label}: copied {copied}, skipped (no comp_details) {skipped}")

print("Done.")
