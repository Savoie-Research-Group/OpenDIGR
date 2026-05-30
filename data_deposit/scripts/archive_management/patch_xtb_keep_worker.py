#!/usr/bin/env python3
"""
Patch existing doi_zips_slim/ zips to add back the two xTB-scan items the
user wants kept:
  - xTB-scan/wbo_timeseries.csv
  - xTB-scan/bem_snapshots/*.json

Reads the source zip in --src-dir (READ-ONLY), pulls just those entries,
writes a new merged zip to <slim>.tmp, then atomic-renames over the slim
zip. Idempotent: skips if the slim already contains both items.

Skips silently if the slim zip does not exist (those are the 118 source
zips that errored during the initial trim and were never produced).
"""
import os, sys, time, argparse, zipfile
from multiprocessing import Pool


def items_to_patch(src_zip):
    """Return {name: bytes} of the entries we want kept from the source."""
    out = {}
    with zipfile.ZipFile(src_zip) as zin:
        for info in zin.infolist():
            n = info.filename
            if n.endswith('/'):
                continue
            if n.endswith('/xTB-scan/wbo_timeseries.csv'):
                out[n] = zin.read(n)
            elif '/xTB-scan/bem_snapshots/' in n:
                out[n] = zin.read(n)
    return out


def patch_one(args):
    src, slim = args
    try:
        # If the slim zip is missing, leave alone (one of the 118 errored zips).
        if not os.path.exists(slim):
            return ("missing_slim", slim, 0)
        if not os.path.exists(src):
            return ("missing_src", src, 0)

        # Read additions from source
        additions = items_to_patch(src)
        if not additions:
            # Source has neither item to add (legitimate; older pipeline)
            return ("no_items", slim, 0)

        # Skip if both items are already present in slim
        with zipfile.ZipFile(slim) as zin:
            existing = set(zin.namelist())
        new = {n: data for n, data in additions.items() if n not in existing}
        if not new:
            return ("already_patched", slim, 0)

        # Build new zip via tmp + atomic rename
        tmp = slim + ".tmp"
        with zipfile.ZipFile(slim) as zin, \
             zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zout:
            for info in zin.infolist():
                zout.writestr(info.filename, zin.read(info.filename))
            for n, data in new.items():
                zout.writestr(n, data)
        os.replace(tmp, slim)
        return ("ok", slim, len(new))
    except Exception as e:
        return ("err", slim, f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keys", required=True)
    ap.add_argument("--src-dir", required=True)
    ap.add_argument("--dst-dir", required=True)
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()

    src_root = os.path.abspath(args.src_dir)
    dst_root = os.path.abspath(args.dst_dir)
    with open(args.keys) as f:
        keys = [ln.strip() for ln in f if ln.strip()]

    pairs = []
    for k in keys:
        suffix = k if k.endswith(".zip") else (k + ".zip")
        pairs.append((os.path.join(src_root, suffix),
                      os.path.join(dst_root, suffix)))

    print(f"Patching {len(pairs)} zips.  src={src_root}  dst={dst_root}  j={args.j}")
    t0 = time.time()
    counts = {"ok":0,"already_patched":0,"missing_slim":0,"missing_src":0,
              "no_items":0,"err":0}
    err_examples = []
    with Pool(args.j) as pool:
        for i,(status,path,info) in enumerate(
                pool.imap_unordered(patch_one, pairs, chunksize=32), 1):
            counts[status] = counts.get(status,0) + 1
            if status == "err" and len(err_examples) < 10:
                err_examples.append((path, info))
            if i % 1000 == 0:
                elapsed = max(time.time()-t0, 0.001)
                rate = i/elapsed
                eta = (len(pairs)-i)/rate if rate else 0
                print(f"  [{i}/{len(pairs)}] {counts}  "
                      f"({rate:.0f}/s, eta {eta:.0f}s)")

    elapsed = time.time()-t0
    print(f"\nDone in {elapsed:.0f}s.  {counts}")
    if err_examples:
        print("Example errors:")
        for p,e in err_examples:
            print(f"  {p}: {e}")


if __name__ == "__main__":
    main()
