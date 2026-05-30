#!/usr/bin/env python3
"""
Trim recipe v2 (2026-05-25, user-revised):
  KEEP: _mat.txt (Mayer/Wiberg/Fuzzy bond-order matrices), .chelpg.xyz
  KEEP: enough of .out to verify calc (input, geometry, SCF, FINAL SP, normal-term)
  DROP: .property.txt, _log.out, xtb-opt/, xTB-scan/
  RECOMPRESS: every entry is deflated (fixes ZIP_STORED injection bug)
"""
import os, re, sys, zipfile

# .out lines we keep (with a generous local window)
KEEP = [
    re.compile(r'^FINAL SINGLE POINT ENERGY'),
    re.compile(r'\*\*\*\*ORCA TERMINATED'),
    re.compile(r'^Total SCF time'),
    re.compile(r'NORMAL TERMINATION'),
    re.compile(r'^! '),                                # input keyword echo
    re.compile(r'^\s*Number of atoms\s+\.+'),
    re.compile(r'^\s*Total Charge\s+'),
    re.compile(r'^\s*Multiplicity\s+'),
    re.compile(r'^CARTESIAN COORDINATES \(ANGSTROEM\)'),
    re.compile(r'^SCF CONVERGENCE'),                   # final SCF block
    re.compile(r'TOTAL RUN TIME'),
]
def slim_out(text):
    L = text.splitlines(); n = len(L); k = [False]*n
    for i, ln in enumerate(L):
        if any(p.search(ln) for p in KEEP):
            for j in range(max(0,i-2), min(n,i+120)): k[j] = True
    return "\n".join(l for l,m in zip(L,k) if m) + "\n"

LOG = re.compile(r'_(mayer|wiberg|fuzzy)_log\.out$')

def trim(src, dst):
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(dst,'w',zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            n = info.filename
            if n.endswith('/'): continue
            if '/xtb-opt/' in n or '/xTB-scan/' in n: continue
            base = os.path.basename(n)
            in_dft = '/DFT-SinglePoint/' in n
            if in_dft:
                # drop .property.txt and Multiwfn logs only
                if base.endswith('.property.txt'): continue
                if LOG.search(base): continue
                data = zin.read(n)
                # slim main DFT .out (not _log.out — those were dropped above)
                if base.endswith('.out'):
                    data = slim_out(data.decode('utf-8','replace')).encode('utf-8')
                zout.writestr(n, data)
            else:
                zout.writestr(n, zin.read(n))

if __name__ == "__main__":
    trim(sys.argv[1], sys.argv[2])
