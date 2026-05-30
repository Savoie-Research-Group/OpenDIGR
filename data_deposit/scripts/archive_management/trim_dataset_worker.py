#!/usr/bin/env python3
"""
Per-chunk worker for the dataset trim job (2026-05-25).

Reads a list of zip keys (paths relative to --src-dir), and for each one:
  - computes the destination at --dst-dir/<key>.zip
  - skips if destination already exists (idempotent re-runs)
  - mkdir -p the destination directory
  - runs the trim recipe from trim_one_zip_v2.py
  - writes to <dst>.tmp then atomic-renames to <dst>.zip so a killed
    task cannot leave a half-written archive in place
  - logs OK/SKIP/ERR per file

Read-only on --src-dir; nothing is ever written there.
"""
import os, sys, time, argparse, importlib.util
from multiprocessing import Pool

from pathlib import Path as _Path
TRIM_SRC = str(_Path(__file__).resolve().parent / "trim_one_zip_v2.py")


def load_trim():
    spec = importlib.util.spec_from_file_location("trim_v2", TRIM_SRC)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.trim


_trim_fn = None
def init_worker():
    global _trim_fn
    _trim_fn = load_trim()


def process_one(args):
    src_path, dst_path = args
    try:
        if os.path.exists(dst_path) and os.path.getsize(dst_path) > 100:
            return ("skip", dst_path, 0)
        os.makedirs(os.path.dirname(dst_path), exist_ok=True)
        tmp = dst_path + ".tmp"
        _trim_fn(src_path, tmp)
        os.replace(tmp, dst_path)
        return ("ok", dst_path, os.path.getsize(dst_path))
    except Exception as e:
        return ("err", src_path, f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keys", required=True)
    ap.add_argument("--src-dir", required=True)
    ap.add_argument("--dst-dir", required=True)
    ap.add_argument("-j", type=int, default=1)
    args = ap.parse_args()

    with open(args.keys) as f:
        keys = [ln.strip() for ln in f if ln.strip()]

    src = os.path.abspath(args.src_dir)
    dst = os.path.abspath(args.dst_dir)

    pairs = []
    for k in keys:
        suffix = k if k.endswith(".zip") else (k + ".zip")
        pairs.append((os.path.join(src, suffix), os.path.join(dst, suffix)))

    print(f"Loaded {len(pairs)} zips. src={src}  dst={dst}  j={args.j}")
    t0 = time.time()
    ok = skip = err = 0
    err_examples = []

    with Pool(args.j, initializer=init_worker) as pool:
        for i, (status, path, info) in enumerate(
                pool.imap_unordered(process_one, pairs, chunksize=32), 1):
            if status == "ok":
                ok += 1
            elif status == "skip":
                skip += 1
            else:
                err += 1
                if len(err_examples) < 10:
                    err_examples.append((path, info))
            if i % 1000 == 0:
                elapsed = max(time.time() - t0, 0.001)
                rate = i / elapsed
                eta = (len(pairs) - i) / rate if rate else 0
                print(f"  [{i}/{len(pairs)}] ok={ok} skip={skip} err={err} "
                      f"({rate:.0f}/s, eta {eta:.0f}s)")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s.  ok={ok}  skip={skip}  err={err}")
    if err_examples:
        print("Example errors:")
        for p, e in err_examples:
            print(f"  {p}: {e}")


if __name__ == "__main__":
    main()
