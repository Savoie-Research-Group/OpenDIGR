#!/usr/bin/env python3
"""Extract functional / basis_set / software / etc. from comp_detail.json files.

Walks doi_files/<DOI>/comp_details/*-comp_detail.json under
GoldDIGR-Comp-details/ and produces:

  results/raw_<field>_per_entry.csv     (one row per JSON file)
  results/raw_<field>_per_paper.csv     (one row per DOI; each distinct
                                         value contributes one vote)
  results/canon_<field>_per_entry.csv   (after canonicalization)
  results/canon_<field>_per_paper.csv
  results/functional_x_basis_per_paper.csv
  results/extraction_log.txt            (null/empty/list/multi-valued counts)
  results/sample_unique_<field>.txt     (first 200 unique raw strings)

Run:
  python3 extract_methods.py
"""

import json
import re
import os
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT      = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/doi_files")
RESULTS   = Path("/groups/bsavoie2/zli43/GoldDIGR-Comp-details/results")
RESULTS.mkdir(exist_ok=True)


# ──────────────────────────────────────────────────────────────────
# Canonicalization rules
# ──────────────────────────────────────────────────────────────────
# Mappings are applied AFTER light normalization (lowercase, strip,
# punctuation collapse). Keys are post-normalization strings.

# light normalization
_PUNCT_RE = re.compile(r"[\s\-_/]+")

def _light(s):
    """Lowercase, strip, collapse whitespace/hyphens/underscores/slashes."""
    if not isinstance(s, str):
        return ""
    s = s.strip().lower()
    s = _PUNCT_RE.sub("", s)
    # remove parentheses contents like "(d,p)" temporarily? no — keep them
    return s

# Canonical-form lookup. The KEY is the post-_light() string.
# The VALUE is the human-readable canonical form we want to display.
# Populate iteratively after inspecting top-100 raw strings.
FUNCTIONAL_CANON = {
    # B3LYP family
    "b3lyp": "B3LYP",
    "ub3lyp": "B3LYP",           # unrestricted variant
    "rb3lyp": "B3LYP",
    "becke3lyp": "B3LYP",
    "becke3lyp3": "B3LYP",
    "b3lypgd2": "B3LYP-D2",
    "b3lypgd3": "B3LYP-D3",
    "b3lypgd3bj": "B3LYP-D3(BJ)",
    "b3lypd3": "B3LYP-D3",
    "b3lypd3bj": "B3LYP-D3(BJ)",
    "b3lypd2": "B3LYP-D2",
    "b3lypd": "B3LYP-D",
    "b3p86": "B3P86",
    "b3pw91": "B3PW91",
    "ub3pw91": "B3PW91",
    # Minnesota functionals
    "m06": "M06",
    "m062x": "M06-2X",
    "m06l": "M06-L",
    "m062xd3": "M06-2X-D3",
    "m06d3": "M06-D3",
    "m05": "M05",
    "m052x": "M05-2X",
    "m11": "M11",
    "m11l": "M11-L",
    "mn15": "MN15",
    "mn15l": "MN15-L",
    "m06hf": "M06-HF",
    # PBE family
    "pbe": "PBE",
    "ggapbe": "PBE",
    "perdewburkeernzerhof": "PBE",
    "perdewburkeernzerhoffunctional": "PBE",
    "pbe1pbe": "PBE0",          # Gaussian's name for PBE0
    "pbe0": "PBE0",
    "pbe0d3": "PBE0-D3",
    "pbed3": "PBE-D3",
    "pbed3bj": "PBE-D3(BJ)",
    "pbed2": "PBE-D2",
    "pbesol": "PBEsol",
    "revpbe": "revPBE",
    # ωB97X family
    "ωb97x": "ωB97X",
    "wb97x": "ωB97X",
    "omegab97x": "ωB97X",
    "ωb97xd": "ωB97X-D",
    "wb97xd": "ωB97X-D",
    "omegab97xd": "ωB97X-D",
    "ωb97xd3": "ωB97X-D3",
    "wb97xd3": "ωB97X-D3",
    "ωb97xv": "ωB97X-V",
    "wb97xv": "ωB97X-V",
    "ωb97m": "ωB97M",
    "wb97m": "ωB97M",
    "ωb97mv": "ωB97M-V",
    "wb97mv": "ωB97M-V",
    # BP86 family
    "bp86": "BP86",
    "bp86d3": "BP86-D3",
    "blyp": "BLYP",
    "blypd3": "BLYP-D3",
    "bpe": "BPE",
    "bp": "BP",
    # TPSS family
    "tpss": "TPSS",
    "tpssh": "TPSSh",
    "tpssd3": "TPSS-D3",
    # CAM-B3LYP
    "camb3lyp": "CAM-B3LYP",
    "camb3lypd3": "CAM-B3LYP-D3",
    # BHandHLYP / HF
    "bhandhlyp": "BHandHLYP",
    "bhandh": "BHandH",
    "hf": "HF",
    "rhf": "HF",
    "uhf": "HF",
    # MN / mPW
    "mpw1pw91": "mPW1PW91",
    "mpw2plyp": "mPW2PLYP",
    "mpwb1k": "mPWB1K",
    # B97 / B97-D
    "b97": "B97",
    "b97d": "B97-D",
    "b97d3": "B97-D3",
    "b97d3bj": "B97-D3(BJ)",
    # Double-hybrid
    "b2plyp": "B2PLYP",
    "b2plypd3": "B2PLYP-D3",
    "dsdpbep86": "DSD-PBEP86",
    # Wave-function
    "ccsdt": "CCSD(T)",
    "ccsd": "CCSD",
    "rccsd": "CCSD",
    "uccsd": "CCSD",
    "mp2": "MP2",
    "rmp2": "MP2",
    "ump2": "MP2",
    "casscf": "CASSCF",
    "caspt2": "CASPT2",
    "nevpt2": "NEVPT2",
    "dlpnoccsdt": "DLPNO-CCSD(T)",
    # X3LYP / LC-ωPBE
    "x3lyp": "X3LYP",
    "lcωpbe": "LC-ωPBE",
    "lcwpbe": "LC-ωPBE",
    # SCAN family
    "scan": "SCAN",
    "r2scan": "r²SCAN",
    "scan0": "SCAN0",
    # HSE
    "hse06": "HSE06",
    "hse": "HSE",
    # PW91
    "pw91": "PW91",
    # RPBE
    "rpbe": "RPBE",
    # GGA
    "gga": "GGA",
    "lda": "LDA",
    # SVWN / LSDA
    "svwn": "SVWN",
    "lsda": "LSDA",
}

BASIS_CANON = {
    # Pople (handle 6-31G* == 6-31G(d), 6-31G** == 6-31G(d,p))
    "631g": "6-31G",
    "631gd": "6-31G(d)",
    "631g*": "6-31G(d)",
    "631gstar": "6-31G(d)",
    "631gdp": "6-31G(d,p)",
    "631g**": "6-31G(d,p)",
    "631gstarstar": "6-31G(d,p)",
    "6311g": "6-311G",
    "6311gd": "6-311G(d)",
    "6311g*": "6-311G(d)",
    "6311gdp": "6-311G(d,p)",
    "6311g**": "6-311G(d,p)",
    "6311gd2p": "6-311G(d,2p)",
    "6311g2df2pd": "6-311G(2df,2pd)",
    "6311gpdp": "6-311+G(d,p)",
    "6311pgdp": "6-311+G(d,p)",
    "6311pg**": "6-311+G(d,p)",
    "6311pg(dp)": "6-311+G(d,p)",
    "6311ppgdp": "6-311++G(d,p)",
    "6311ppg**": "6-311++G(d,p)",
    "6311ppg(dp)": "6-311++G(d,p)",
    "6311ppg2dp": "6-311++G(2d,p)",
    "6311ppg2d2p": "6-311++G(2d,2p)",
    "6311pg2dp": "6-311+G(2d,p)",
    "6311pg2dp:": "6-311+G(2d,p)",
    "6311pg(2dp)": "6-311+G(2d,p)",
    "6311ppg(dpd)": "6-311++G(d,p)",
    "631ppgdp": "6-31++G(d,p)",
    "631pgd": "6-31+G(d)",
    "631pg*": "6-31+G(d)",
    "631pgdp": "6-31+G(d,p)",
    "631pg**": "6-31+G(d,p)",
    "321g": "3-21G",
    "321g*": "3-21G(d)",
    "sto3g": "STO-3G",
    # Effective Core Potentials
    "lanl2dz": "LANL2DZ",
    "lanl2tz": "LANL2TZ",
    "lanl08": "LANL08",
    "sdd": "SDD",
    "stuttgart": "Stuttgart",
    "stuttgartrsc1997": "Stuttgart RSC 1997",
    # Karlsruhe (def2)
    "def2sv": "def2-SV",
    "def2svp": "def2-SVP",
    "def2svpd": "def2-SVPD",
    "def2tzvp": "def2-TZVP",
    "def2tzvpp": "def2-TZVPP",
    "def2tzvpd": "def2-TZVPD",
    "def2qzvp": "def2-QZVP",
    "def2qzvpp": "def2-QZVPP",
    # Dunning
    "ccpvdz": "cc-pVDZ",
    "ccpvtz": "cc-pVTZ",
    "ccpvqz": "cc-pVQZ",
    "ccpv5z": "cc-pV5Z",
    "augccpvdz": "aug-cc-pVDZ",
    "augccpvtz": "aug-cc-pVTZ",
    "augccpvqz": "aug-cc-pVQZ",
    # Plane-wave (all variants)
    "planewavebasisset": "plane-wave",
    "planewave": "plane-wave",
    "planewaves": "plane-wave",
    "pw": "plane-wave",
    "paw": "PAW",
    # Ahlrichs
    "tzv": "TZV",
    "tzvp": "TZVP",
    "tzvpp": "TZVPP",
    "tz2p": "TZ2P",
    "dz": "DZ",
    "dzp": "DZP",
    "qzvp": "QZVP",
    # ma- (minimally augmented Karlsruhe)
    "ma6311gdp": "ma-6-311G(d,p)",
    "ma6311ppgdp": "ma-6-311++G(d,p)",
    # mixed / unspecified
    "mixed": "(mixed basis)",
    "mixedbasis": "(mixed basis)",
    "mixedbasisset": "(mixed basis)",
}

SOFTWARE_CANON = {
    "gaussian": "Gaussian",
    "gaussian03": "Gaussian 03",
    "gaussian09": "Gaussian 09",
    "gaussian09w": "Gaussian 09",
    "gaussian16": "Gaussian 16",
    "gaussian98": "Gaussian 98",
    "g03": "Gaussian 03",
    "g09": "Gaussian 09",
    "g16": "Gaussian 16",
    "g98": "Gaussian 98",
    "orca": "ORCA",
    "vasp": "VASP",
    "viennaabinitiosimulationpackage": "VASP",
    "viennaabinitiosimulationpackagevasp": "VASP",
    "qchem": "Q-Chem",
    "molpro": "Molpro",
    "turbomole": "Turbomole",
    "nwchem": "NWChem",
    "psi4": "Psi4",
    "amber": "AMBER",
    "gromacs": "GROMACS",
    "lammps": "LAMMPS",
    "cp2k": "CP2K",
    "siesta": "SIESTA",
    "adf": "ADF",
    "amsadf": "AMS/ADF",
    "ams": "AMS",
    "dalton": "Dalton",
    "jaguar": "Jaguar",
    "spartan": "Spartan",
    "crystal": "CRYSTAL",
    "crystal17": "CRYSTAL17",
    "crystal14": "CRYSTAL14",
    "quantumespresso": "Quantum ESPRESSO",
    "quantumespresso(qe)": "Quantum ESPRESSO",
    "qe": "Quantum ESPRESSO",
    "castep": "CASTEP",
    "multiwfn": "Multiwfn",
    "materialsstudio": "Materials Studio",
    "dmol3": "DMol3",
    "vmd": "VMD",
    "namd": "NAMD",
    "abinit": "ABINIT",
    "wien2k": "WIEN2k",
    "siesta4": "SIESTA",
    "fhiaims": "FHI-aims",
    "openmolcas": "OpenMolcas",
    "molcas": "Molcas",
    "ezdoublec": "EzdoubleC",
    "pyscf": "PySCF",
    "cfour": "CFOUR",
    "deepmd": "DeePMD",
    "lobster": "LOBSTER",
}


# Software version-stripping regex.  If a software string starts with a
# known root (Gaussian, ORCA, VASP, etc.) followed by a version/revision
# token, collapse it to the root.  Applied AFTER canonicalize lookup.
_SOFTWARE_ROOTS = [
    ("Gaussian 16", re.compile(r"^Gaussian\s*16\b", re.IGNORECASE)),
    ("Gaussian 09", re.compile(r"^Gaussian\s*0?9\b", re.IGNORECASE)),
    ("Gaussian 03", re.compile(r"^Gaussian\s*0?3\b", re.IGNORECASE)),
    ("Gaussian 98", re.compile(r"^Gaussian\s*98\b", re.IGNORECASE)),
    ("Gaussian",    re.compile(r"^Gaussian\b", re.IGNORECASE)),
    ("ORCA",        re.compile(r"^ORCA\b", re.IGNORECASE)),
    ("VASP",        re.compile(r"^VASP\b", re.IGNORECASE)),
    ("VASP",        re.compile(r"^Vienna Ab[\s-]?initio", re.IGNORECASE)),
    ("CASTEP",      re.compile(r"^CASTEP\b", re.IGNORECASE)),
    ("Turbomole",   re.compile(r"^Turbomole\b", re.IGNORECASE)),
    ("Q-Chem",      re.compile(r"^Q[-\s]?Chem\b", re.IGNORECASE)),
    ("NWChem",      re.compile(r"^NWChem\b", re.IGNORECASE)),
    ("Molpro",      re.compile(r"^Molpro\b", re.IGNORECASE)),
    ("CP2K",        re.compile(r"^CP2K\b", re.IGNORECASE)),
    ("Quantum ESPRESSO", re.compile(r"^Quantum\s*ESPRESSO\b", re.IGNORECASE)),
    ("Quantum ESPRESSO", re.compile(r"^QE\b", re.IGNORECASE)),
    ("ADF",         re.compile(r"^ADF\b", re.IGNORECASE)),
    ("AMS",         re.compile(r"^AMS\b", re.IGNORECASE)),
    ("DMol3",       re.compile(r"^DMol3?\b", re.IGNORECASE)),
    ("Jaguar",      re.compile(r"^Jaguar\b", re.IGNORECASE)),
    ("Materials Studio", re.compile(r"^Materials\s*Studio", re.IGNORECASE)),
    ("CRYSTAL",     re.compile(r"^CRYSTAL\d*\b", re.IGNORECASE)),
    ("Spartan",     re.compile(r"^Spartan\b", re.IGNORECASE)),
    ("WIEN2k",      re.compile(r"^WIEN2k\b", re.IGNORECASE)),
    ("FHI-aims",    re.compile(r"^FHI[-\s]?aims\b", re.IGNORECASE)),
    ("OpenMolcas",  re.compile(r"^OpenMolcas\b", re.IGNORECASE)),
    ("Molcas",      re.compile(r"^Molcas\b", re.IGNORECASE)),
    ("Multiwfn",    re.compile(r"^Multiwfn\b", re.IGNORECASE)),
]

def collapse_software(s):
    """Map e.g. 'Gaussian 09 Revision D.01' -> 'Gaussian 09', 'ORCA 5.0.3' -> 'ORCA'."""
    if not isinstance(s, str):
        return s
    for canonical, pattern in _SOFTWARE_ROOTS:
        if pattern.search(s.strip()):
            return canonical
    return s


def canonicalize(raw, table):
    """raw -> canonical via lookup, fall back to raw if not in table."""
    key = _light(raw)
    return table.get(key, raw.strip() if isinstance(raw, str) else "")


def canonicalize_software(raw):
    """Two-step: table lookup, then version-stripping regex collapse."""
    s = canonicalize(raw, SOFTWARE_CANON)
    return collapse_software(s)


# ──────────────────────────────────────────────────────────────────
# Walking + extraction
# ──────────────────────────────────────────────────────────────────

def doi_from_path(p):
    """Pull '10.xxxx/yyyy' from a path like
    .../doi_files/10.xxxx/yyyy/comp_details/<file>.json"""
    parts = p.parts
    for i, part in enumerate(parts):
        if part.startswith("10.") and i + 1 < len(parts):
            return part + "/" + parts[i + 1]
    return "unknown"


def normalize_value(v):
    """Convert a field value into a list of non-empty strings."""
    if v is None:
        return []
    if isinstance(v, list):
        out = []
        for item in v:
            if isinstance(item, str) and item.strip():
                out.append(item.strip())
        return out
    if isinstance(v, str):
        v = v.strip()
        return [v] if v else []
    return []


def main():
    json_files = list(ROOT.rglob("*-comp_detail.json"))
    n_total = len(json_files)
    print("Found {} comp_detail.json files".format(n_total), file=sys.stderr)

    # Counters
    fields = ["software", "basis_set", "functional",
              "level_of_theory", "special_treatments",
              "method_family", "calculation_types"]

    raw_entry_counts  = {f: Counter() for f in fields}
    raw_paper_seen    = {f: defaultdict(set) for f in fields}  # doi -> set of values

    # Special: functional × basis combo (per paper)
    combo_paper_seen  = defaultdict(set)  # doi -> set of (func, basis)

    # Log counters
    log = Counter()
    parse_errors = []

    for i, jp in enumerate(json_files):
        if i % 5000 == 0:
            print("  {}/{}".format(i, n_total), file=sys.stderr)
        try:
            with open(jp) as f:
                obj = json.load(f)
        except Exception as e:
            log["parse_error"] += 1
            parse_errors.append((str(jp), str(e)))
            continue

        # A few JSONs are a list of dicts (multiple calculations per SI).
        # Treat each dict as an independent entry under the same DOI.
        if isinstance(obj, list):
            log["list_typed_file"] += 1
            entries = [o for o in obj if isinstance(o, dict)]
        elif isinstance(obj, dict):
            entries = [obj]
        else:
            log["unexpected_type"] += 1
            continue

        doi = doi_from_path(jp)

        for obj in entries:
            log["entries"] += 1

            # Per-paper functional × basis
            funcs   = normalize_value(obj.get("functional"))
            bases   = normalize_value(obj.get("basis_set"))
            for fn in funcs:
                for bs in bases:
                    cf = canonicalize(fn, FUNCTIONAL_CANON)
                    cb = canonicalize(bs, BASIS_CANON)
                    combo_paper_seen[doi].add((cf, cb))

            for field in fields:
                vals = normalize_value(obj.get(field))
                if not vals:
                    log["empty_" + field] += 1
                    continue
                for v in vals:
                    raw_entry_counts[field][v] += 1
                    raw_paper_seen[field][doi].add(v)

    # ── per-paper raw counts: each DOI contributes one vote per distinct value
    raw_paper_counts = {f: Counter() for f in fields}
    for field in fields:
        for doi, vals in raw_paper_seen[field].items():
            for v in vals:
                raw_paper_counts[field][v] += 1

    # ── canonical counts
    canon_entry_counts = {f: Counter() for f in fields}
    canon_paper_counts = {f: Counter() for f in fields}
    canon_tables = {
        "functional":         FUNCTIONAL_CANON,
        "basis_set":          BASIS_CANON,
        "software":           SOFTWARE_CANON,
        "level_of_theory":    {},
        "special_treatments": {},
        "method_family":      {},
        "calculation_types":  {},
    }
    for field in fields:
        table = canon_tables[field]
        # software gets the extra version-stripping pass
        if field == "software":
            canon_fn = canonicalize_software
        else:
            canon_fn = lambda v, t=table: canonicalize(v, t)
        for v, c in raw_entry_counts[field].items():
            canon_entry_counts[field][canon_fn(v)] += c
        for doi, vals in raw_paper_seen[field].items():
            seen = set(canon_fn(v) for v in vals)
            for cv in seen:
                canon_paper_counts[field][cv] += 1

    # ── write CSVs
    def write_hist(path, counter):
        rows = counter.most_common()
        with open(path, "w", newline="") as fp:
            w = csv.writer(fp)
            w.writerow(["value", "count"])
            for v, c in rows:
                w.writerow([v, c])

    for field in fields:
        write_hist(RESULTS / "raw_{}_per_entry.csv".format(field),
                   raw_entry_counts[field])
        write_hist(RESULTS / "raw_{}_per_paper.csv".format(field),
                   raw_paper_counts[field])
        write_hist(RESULTS / "canon_{}_per_entry.csv".format(field),
                   canon_entry_counts[field])
        write_hist(RESULTS / "canon_{}_per_paper.csv".format(field),
                   canon_paper_counts[field])

    # functional × basis combo (per paper, canonical)
    combo_paper_counts = Counter()
    for doi, combos in combo_paper_seen.items():
        for fn, bs in combos:
            combo_paper_counts[(fn, bs)] += 1
    with open(RESULTS / "functional_x_basis_per_paper.csv", "w", newline="") as fp:
        w = csv.writer(fp)
        w.writerow(["functional", "basis_set", "papers"])
        for (fn, bs), c in combo_paper_counts.most_common():
            w.writerow([fn, bs, c])

    # extraction log
    with open(RESULTS / "extraction_log.txt", "w") as fp:
        fp.write("total comp_detail.json files: {}\n".format(n_total))
        fp.write("successfully parsed: {}\n".format(log["entries"]))
        fp.write("parse errors: {}\n".format(log["parse_error"]))
        fp.write("\nper-field empty/null counts:\n")
        for f in fields:
            fp.write("  {:24s} empty: {}\n".format(f, log["empty_" + f]))
        fp.write("\nunique values per field (raw):\n")
        for f in fields:
            fp.write("  {:24s} {} unique\n".format(f, len(raw_entry_counts[f])))
        fp.write("\nunique papers seen: {}\n".format(
            len(set().union(*[set(raw_paper_seen[f].keys()) for f in fields]))))
        if parse_errors:
            fp.write("\nfirst 10 parse errors:\n")
            for p, e in parse_errors[:10]:
                fp.write("  {}: {}\n".format(p, e))

    # sample unique strings per field (for canonicalization iteration)
    for field in fields:
        with open(RESULTS / "sample_unique_{}.txt".format(field), "w") as fp:
            for v, c in raw_entry_counts[field].most_common(500):
                fp.write("{:>8d}\t{}\n".format(c, v))

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
