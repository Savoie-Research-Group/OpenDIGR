#!/usr/bin/env python3
"""
regen_frame_jsons.py

Regenerate IRC_Analysis/<subdir>/frame_NNNNN.json files inside a tar.zst
archive using the *patched* yarp. Schema matches the originals:

  {"source": "<subdir>/frame_<NNNNN>",
   "atoms": [{"id": "Ir0", "el": "Ir", "e": <int>}, ...],
   "bonds": [{"i": "Ir0", "j": "C1", "order": <int>}, ...]}

Only re-writes JSONs that are present in the archive. CSVs (which doi_tar_zsts
doesn't store) are not produced -- those live in doi_zips and will be handled
in phase 2.

Charge comes from the filename suffix (e.g. 25_0_1.tar.zst -> charge=0).
yarpecule is called with canon=False so atom indices match the source xyz.

Usage:
  python3 regen_frame_jsons.py <archive.tar.zst>
      Validates and writes <archive>.fixed.tar.zst alongside the input.
  python3 regen_frame_jsons.py <archive.tar.zst> --inplace
      Overwrites the archive after writing to .fixed.tmp and renaming.
  python3 regen_frame_jsons.py <archive.tar.zst> --minimal
      Only regenerate finished_first + finished_last frame_00000.json
      (the two files extract_oxidation_states.py reads). Skip the IRC
      sweep, initial-TS and TS subdirs.
"""
import argparse, contextlib, io, json, os, re, subprocess, sys, tarfile, tempfile
from pathlib import Path

import numpy as np
import yarp as yp

CHARGE_RE = re.compile(r".*_(-?\d+)_(\d+)\.tar\.zst$")

# Map IRC_Analysis subdir -> source file inside the archive
# (finished_irc reads multi-frame trj; others read single-frame xyz).
SUBDIR_SOURCES_FULL = {
    "finished_first": ("finished_first.xyz", "single"),
    "finished_last":  ("finished_last.xyz",  "single"),
    "finished_irc":   ("finished_irc.trj",   "multi"),
    "initial-TS":     ("input.xyz",          "single"),
    "TS":             ("ts_final_geometry.xyz", "single"),
}
SUBDIR_SOURCES_MINIMAL = {
    "finished_first": ("finished_first.xyz", "single"),
    "finished_last":  ("finished_last.xyz",  "single"),
}


def silent_yarpecule(xyz_text: str, charge: int):
    """Run yarpecule on xyz given as text. Writes to temp file, returns yarpecule."""
    with tempfile.NamedTemporaryFile("w", suffix=".xyz", delete=False) as tf:
        tf.write(xyz_text)
        tf.flush()
        tmp_path = tf.name
    try:
        with open(os.devnull, "w") as dn, \
             contextlib.redirect_stdout(dn), contextlib.redirect_stderr(dn):
            y = yp.yarpecule((tmp_path, charge), canon=False)
    finally:
        try: os.unlink(tmp_path)
        except OSError: pass
    return y


def parse_trj(text: str):
    """Yield (frame_index, frame_xyz_text) for each frame in a multi-frame xyz/trj."""
    lines = text.splitlines()
    i = 0
    frame_idx = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        try:
            n = int(lines[i].strip())
        except ValueError:
            break
        # one frame: header(1) + comment(1) + n atoms
        block = lines[i : i + 2 + n]
        if len(block) < 2 + n:
            break
        yield frame_idx, "\n".join(block) + "\n"
        i += 2 + n
        frame_idx += 1


def build_json(yarpecule, source_label: str) -> dict:
    """Build the doi_tar_zsts frame JSON from a yarpecule's first bond_mat."""
    elements = [str(e).capitalize() for e in yarpecule.elements]
    n = len(elements)
    be = yarpecule.bond_mats[0]
    adj = yarpecule.adj_mat
    ids = [f"{el}{i}" for i, el in enumerate(elements)]
    atoms = [{"id": ids[i], "el": elements[i], "e": int(be[i, i])} for i in range(n)]
    bonds = []
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j]:
                bonds.append({"i": ids[i], "j": ids[j], "order": int(be[i, j])})
    return {"source": source_label, "atoms": atoms, "bonds": bonds}


def regenerate_archive(path: Path, output: Path, minimal: bool = False):
    sources = SUBDIR_SOURCES_MINIMAL if minimal else SUBDIR_SOURCES_FULL
    m = CHARGE_RE.match(path.name)
    if m is None:
        raise SystemExit(f"Cannot parse charge from filename: {path.name}")
    charge = int(m.group(1))

    # Decompress entire archive into memory
    proc = subprocess.run(["zstd", "-dc", str(path)], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True)
    raw = proc.stdout

    members = {}
    with tarfile.open(fileobj=io.BytesIO(raw), mode="r:") as tf:
        for ti in tf.getmembers():
            if ti.isreg():
                f = tf.extractfile(ti)
                members[ti.name] = (ti, f.read() if f else b"")
            else:
                members[ti.name] = (ti, b"")

    # Cache: source_path -> { frame_idx: yarpecule_json_dict }
    cache = {}
    n_regen = 0
    n_skipped = 0
    log = []

    for name in list(members):
        # Only process IRC_Analysis/<subdir>/frame_NNNNN.json files
        parts = name.split("/")
        if len(parts) < 3 or parts[0] != "IRC_Analysis":
            continue
        subdir = parts[1]
        if subdir not in sources:
            continue
        fname = parts[2]
        m2 = re.fullmatch(r"frame_(\d{5})\.json", fname)
        if not m2:
            continue
        frame_idx = int(m2.group(1))

        src_name, mode = sources[subdir]
        if src_name not in members:
            log.append(f"  SKIP {name}: source {src_name} not in archive")
            n_skipped += 1
            continue

        # Build/lookup yarpecule for this source+frame
        cache_key = (src_name, frame_idx)
        if cache_key not in cache:
            src_text = members[src_name][1].decode("utf-8", errors="replace")
            if mode == "single":
                if frame_idx != 0:
                    log.append(f"  SKIP {name}: single-frame source but frame_idx={frame_idx}")
                    n_skipped += 1
                    continue
                xyz_text = src_text if src_text.strip().split("\n", 1)[0].strip().isdigit() \
                          else src_text  # already an xyz
                y = silent_yarpecule(xyz_text, charge)
                cache[cache_key] = y
            else:  # multi
                want = None
                for fi, xyz_text in parse_trj(src_text):
                    if fi == frame_idx:
                        want = xyz_text
                        break
                if want is None:
                    log.append(f"  SKIP {name}: frame {frame_idx} not found in {src_name}")
                    n_skipped += 1
                    continue
                y = silent_yarpecule(want, charge)
                cache[cache_key] = y

        y = cache[cache_key]
        source_label = f"{subdir}/frame_{frame_idx:05d}"
        new_json = build_json(y, source_label)
        new_bytes = json.dumps(new_json, separators=(",", ":")).encode("utf-8")

        # Update member in-place (mutate size + buffer)
        ti = members[name][0]
        ti.size = len(new_bytes)
        members[name] = (ti, new_bytes)
        n_regen += 1

    # Repack tar (preserve original order)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:") as tf:
        for ti, data in members.values():
            if ti.isreg():
                tf.addfile(ti, io.BytesIO(data))
            else:
                tf.addfile(ti)

    # Compress with zstd
    tar_bytes = buf.getvalue()
    proc = subprocess.run(["zstd", "-T0", "-19", "-c"], input=tar_bytes,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    output.write_bytes(proc.stdout)

    print(f"  regenerated: {n_regen}  skipped: {n_skipped}", file=sys.stderr)
    for ln in log:
        print(ln, file=sys.stderr)
    return n_regen, n_skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archive", type=Path)
    ap.add_argument("--inplace", action="store_true")
    ap.add_argument("--minimal", action="store_true",
                    help="Only regenerate finished_first + finished_last frame_00000.json.")
    args = ap.parse_args()
    if args.inplace:
        tmp = args.archive.with_suffix(args.archive.suffix + ".tmp")
        regenerate_archive(args.archive, tmp, minimal=args.minimal)
        tmp.replace(args.archive)
        print(f"Wrote in place: {args.archive}", file=sys.stderr)
    else:
        out = args.archive.with_suffix(".fixed.tar.zst")
        if args.archive.name.endswith(".tar.zst"):
            out = args.archive.parent / (args.archive.name[:-len(".tar.zst")] + ".fixed.tar.zst")
        regenerate_archive(args.archive, out, minimal=args.minimal)
        print(f"Wrote: {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
