#!/usr/bin/env python3
"""
regen_zip_archive.py <path/to/archive.zip>

Full in-place regeneration of one doi_zips_slim archive:
  1) unzip to /tmp/regen_<pid>/<reaction>/
  2) delete the existing yarp-derived files (frame_*.{json,csv} and
     reaction/*.clean.json and sankey/*); xtb-scan, DFT-SinglePoint,
     molsimp are NOT touched
  3) run yarp_results_builder.process_directory(overwrite=True), which
     rewrites all frame_*.json + _bond_electrons.csv + _adjacency.csv files
     AND reaction/deterministic_{forward,reverse}.clean.json
  4) run irc_sankey.py --root <IRC_Analysis/finished_irc> --out <IRC_Analysis/sankey>
  5) repack zip atomically (<archive>.tmp -> rename)

The patched YARP is picked up automatically via the env's installed yarp.
"""
from __future__ import annotations
import argparse, os, re, shutil, subprocess, sys, tempfile, time, zipfile
from pathlib import Path

# Add the metal_ligand_pkg_clean directory to PYTHONPATH so the import works
METAL_LIGAND_PKG = "/home/li1724/metal_ligand_pkg_clean"
sys.path.insert(0, METAL_LIGAND_PKG)

# Import yarp_results_builder's main worker
import yarp_results_builder as yrb

SANKEY_SCRIPT = os.path.join(METAL_LIGAND_PKG, "irc_sankey.py")

# Filename-encoded charge: e.g. 25_0_1.zip -> charge=0; TS-Frag-CyXantCHD-ax_-1_1.zip -> charge=-1
CHARGE_RE = re.compile(r".*_(-?\d+)_(\d+)$")


def parse_charge(stem: str) -> int | None:
    m = CHARGE_RE.match(stem)
    if m is None:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def regen_one(zip_path: Path, keep_temp: bool = False, run_sankey: bool = True,
              verbose: bool = True) -> tuple[bool, str]:
    t0 = time.time()
    work = Path(tempfile.mkdtemp(prefix=f"regen_{os.getpid()}_", dir="/tmp"))
    try:
        # 1) Unzip
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(work)

        # The zip should contain a single top-level dir matching the reaction stem
        entries = [p for p in work.iterdir() if p.is_dir()]
        if len(entries) != 1:
            return False, f"unexpected entries in zip: {[p.name for p in entries]}"
        reaction_dir = entries[0]
        stem = reaction_dir.name
        charge = parse_charge(stem)

        # 2) Wipe yarp-derived files inside IRC_Analysis (keep DFT-SinglePoint,
        #    xTB-scan, molsimp, top-level xyz/trj/yaml/etc untouched)
        irc_root = reaction_dir / "IRC_Analysis"
        if irc_root.exists():
            for sub in ("finished_first","finished_last","finished_irc","initial-TS","TS"):
                d = irc_root / sub
                if d.exists():
                    for f in d.iterdir():
                        # keep ts_frame.txt — yarp_results_builder will rewrite it
                        if f.is_file():
                            f.unlink()
            # Wipe only deterministic_{forward,reverse}.clean.json — leave
            # Fail*-deterministic_*.clean.json (those are produced by separate
            # FailXX checker scripts and are out of scope for this OS regen).
            reaction_d = irc_root / "reaction"
            if reaction_d.exists():
                for fname in ("deterministic_forward.clean.json",
                              "deterministic_reverse.clean.json"):
                    p = reaction_d / fname
                    if p.exists():
                        p.unlink()
            # Wipe sankey/ entirely (it will be fully rewritten)
            sankey_d = irc_root / "sankey"
            if sankey_d.exists():
                for f in sankey_d.iterdir():
                    if f.is_file():
                        f.unlink()

        # 3) Run yarp_results_builder.process_directory (overwrite=True)
        # GUARD (2026-05-31 ZL): yarp's table_generator calls quit() (SystemExit)
        # when input.xyz has unknown elements (e.g. uppercase 'TI'). yrb only
        # catches Exception per frame, so SystemExit propagates and kills the
        # script BEFORE the rezip step → the in-place archive silently retains
        # its pre-patch JSONs. Catch BaseException so we still rezip whatever
        # finished_first/last/irc frames yrb did write before dying.
        yrb_warn = ""
        try:
            ok = yrb.process_directory(reaction_dir,
                                       overwrite=True,
                                       verbose=verbose,
                                       charge=charge)
        except BaseException as e:
            ok = True
            yrb_warn = f" yrb_aborted={type(e).__name__}"
            if verbose:
                print(f"  -> WARN: yrb raised {type(e).__name__}: {e} — continuing to rezip with partial regen",
                      file=sys.stderr)
        if not ok:
            return False, "yarp_results_builder.process_directory returned False"

        # 4) Run irc_sankey.py — note: --root must be IRC_Analysis (it appends
        #    /finished_irc internally to find ts_frame.txt and per-frame data).
        if run_sankey:
            finished_irc_dir = irc_root / "finished_irc"
            sankey_out = irc_root / "sankey"
            if finished_irc_dir.exists() and (finished_irc_dir / "ts_frame.txt").exists():
                sankey_out.mkdir(parents=True, exist_ok=True)
                proc = subprocess.run(
                    [sys.executable, SANKEY_SCRIPT,
                     "--root", str(irc_root),
                     "--out", str(sankey_out),
                     "--window", "20"],
                    capture_output=True, text=True, timeout=120)
                if proc.returncode != 0 and verbose:
                    print(f"  -> warning: irc_sankey rc={proc.returncode}", file=sys.stderr)
                    print(proc.stderr[-500:], file=sys.stderr)

        # 5) Repack zip atomically — write to <archive>.tmp, then rename
        tmp_zip = zip_path.with_suffix(zip_path.suffix + ".tmp")
        if tmp_zip.exists():
            tmp_zip.unlink()
        with zipfile.ZipFile(tmp_zip, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zo:
            for root, dirs, files in os.walk(reaction_dir):
                root_p = Path(root)
                for fname in files:
                    fp = root_p / fname
                    # arcname rooted at the reaction folder name
                    arcname = fp.relative_to(work)
                    zo.write(fp, arcname=str(arcname))
        os.replace(tmp_zip, zip_path)

        return True, f"ok in {time.time()-t0:.1f}s{yrb_warn}"
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"
    finally:
        if not keep_temp:
            shutil.rmtree(work, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archive", type=Path)
    ap.add_argument("--no-sankey", action="store_true")
    ap.add_argument("--keep-temp", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    ok, msg = regen_one(args.archive,
                        keep_temp=args.keep_temp,
                        run_sankey=not args.no_sankey,
                        verbose=not args.quiet)
    tag = "OK" if ok else "FAIL"
    print(f"{tag}  {args.archive}  {msg}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
