#!/usr/bin/env python3
"""
detect_cc_metadata.py  –  flag text files that contain computational‑chemistry metadata

Usage
-----
    python detect_cc_metadata.py                # scan current directory tree
    python detect_cc_metadata.py /path/to/dir   # scan another directory

Output
------
Creates (or overwrites) a file called  chem_metadata_flagged.txt
containing one absolute path per line for every matching .txt file found.
"""

import re
import sys
from pathlib import Path

# ───── keyword / phrase list ─────────────────────────────────────────────────────
PATTERNS = [
    r"\bcomputational\b",
    r"\btheoretical\b",
    r"\bDFT\b",
    r"\bcoordinates?\b",
    r"\bGaussian\b",
    r"\bORCA\b",
    r"\bJaguar\b",
    r"\boptimization\b",
    r"\bbasis\s+set\b",
    r"\bfunctional\b",
    r"\bCCSD\b",
    r"\bcoupled\s+cluster\b",
    r"\bPBE\b",
    r"\bB3LYP\b",
    r"\b6-31G\b",
    r"\bM06\b",
    r"\bCPCM\b",
    r"\bsolvation\s+effect\b",
    r"\bdispersion\b",
    r"\bpolarization\b",
    r"\bfrequency\b",
    r"\bimaginary\b",
    r"\btransition\s+state\b",
]

# Combine into one case‑insensitive regex (OR’ed with |)
META_REGEX = re.compile("|".join(PATTERNS), re.IGNORECASE)

# ───── helper : does file contain any metadata patterns? ─────────────────────────
def has_metadata(txt_path: Path) -> bool:
    """Return True if any pattern is found in the file’s text."""
    try:
        with txt_path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if META_REGEX.search(line):
                    return True
    except Exception as exc:
        # Could not read the file – silently skip but print a warning
        print(f"[warn] Cannot read {txt_path}: {exc}", file=sys.stderr)
    return False


# ───── main ──────────────────────────────────────────────────────────────────────
def main(root: Path) -> None:
    parent = root.parent
    name = root.name
    Path(f"{parent}/comp_details").mkdir(exist_ok=True)
    flagged_file = Path(f"{parent}/comp_details/{name}-chem_metadata_flagged.txt").resolve()
    if flagged_file.is_file():
        print(f"{flagged_file} is already there. Skip")
        return
    matches = []

    for txt in root.rglob("*.txt"):          # recursive *.txt search
        if txt.is_file() and has_metadata(txt):
            matches.append(txt.resolve())

    # Write (or overwrite) results
    flagged_file.write_text("\n".join(map(str, matches)) + ("\n" if matches else ""))
    print(f"✓ Scanned '{root}'.  {len(matches)} file(s) flagged → {flagged_file}")

if __name__ == "__main__":
    scan_root = Path(sys.argv[1]).expanduser().resolve() if len(sys.argv) > 1 else Path.cwd()
    main(scan_root)

