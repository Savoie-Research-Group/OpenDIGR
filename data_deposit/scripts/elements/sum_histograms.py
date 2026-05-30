#!/usr/bin/env python3
# sum_histograms.py
"""
Sum per-process element histograms and optionally write a PT heatmap SVG.

Input CSV format expected (header present):
  element,count_files_with_element

Usage:
  # Sum all CSVs in a directory:
  python sum_histograms.py --csv-glob "hist_out/elements_hist_*.csv" \
                           --out-csv final_element_counts.csv \
                           --svg final_pt_heatmap.svg

If --svg is omitted, only the summed CSV is produced.
"""

import argparse, csv, glob, math, os
from collections import Counter

# Periodic table layout (group 1..18, period) for main blocks + Ln/An rows.
# Minimal mapping (extend if you care about superheavy placement tweaks):
PT_POS = {
#  group, period
'H':(1,1),'He':(18,1),
'Li':(1,2),'Be':(2,2),'B':(13,2),'C':(14,2),'N':(15,2),'O':(16,2),'F':(17,2),'Ne':(18,2),
'Na':(1,3),'Mg':(2,3),'Al':(13,3),'Si':(14,3),'P':(15,3),'S':(16,3),'Cl':(17,3),'Ar':(18,3),
'K':(1,4),'Ca':(2,4),'Sc':(3,4),'Ti':(4,4),'V':(5,4),'Cr':(6,4),'Mn':(7,4),'Fe':(8,4),'Co':(9,4),'Ni':(10,4),'Cu':(11,4),'Zn':(12,4),
'Ga':(13,4),'Ge':(14,4),'As':(15,4),'Se':(16,4),'Br':(17,4),'Kr':(18,4),
'Rb':(1,5),'Sr':(2,5),'Y':(3,5),'Zr':(4,5),'Nb':(5,5),'Mo':(6,5),'Tc':(7,5),'Ru':(8,5),'Rh':(9,5),'Pd':(10,5),'Ag':(11,5),'Cd':(12,5),
'In':(13,5),'Sn':(14,5),'Sb':(15,5),'Te':(16,5),'I':(17,5),'Xe':(18,5),
'Cs':(1,6),'Ba':(2,6),
'La':(3,6),'Ce':(4,8),'Pr':(5,8),'Nd':(6,8),'Pm':(7,8),'Sm':(8,8),'Eu':(9,8),'Gd':(10,8),'Tb':(11,8),'Dy':(12,8),'Ho':(13,8),'Er':(14,8),'Tm':(15,8),'Yb':(16,8),'Lu':(17,8),
'Hf':(4,6),'Ta':(5,6),'W':(6,6),'Re':(7,6),'Os':(8,6),'Ir':(9,6),'Pt':(10,6),'Au':(11,6),'Hg':(12,6),
'Tl':(13,6),'Pb':(14,6),'Bi':(15,6),'Po':(16,6),'At':(17,6),'Rn':(18,6),
'Fr':(1,7),'Ra':(2,7),
'Ac':(3,7),'Th':(4,9),'Pa':(5,9),'U':(6,9),'Np':(7,9),'Pu':(8,9),'Am':(9,9),'Cm':(10,9),'Bk':(11,9),'Cf':(12,9),'Es':(13,9),'Fm':(14,9),'Md':(15,9),'No':(16,9),'Lr':(17,9),
'Rf':(4,7),'Db':(5,7),'Sg':(6,7),'Bh':(7,7),'Hs':(8,7),'Mt':(9,7),'Ds':(10,7),'Rg':(11,7),'Cn':(12,7),'Nh':(13,7),'Fl':(14,7),'Mc':(15,7),'Lv':(16,7),'Ts':(17,7),'Og':(18,7)
}

# Order for writing CSV nicely
Z_ORDER = [
'H','He',
'Li','Be','B','C','N','O','F','Ne',
'Na','Mg','Al','Si','P','S','Cl','Ar',
'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',
'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'
]

LAST = "Rn"
if LAST in Z_ORDER:
    cutoff = Z_ORDER.index(LAST) + 1
    Z_ORDER = Z_ORDER[:cutoff]
    ALLOWED = set(Z_ORDER)
    PT_POS = {el: pos for el, pos in PT_POS.items() if el in ALLOWED}
else:
    ALLOWED = set(Z_ORDER)  # fallback (shouldn’t happen)

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv-glob", required=True, help='Glob for per-process CSVs, e.g. "hist_out/elements_hist_*.csv".')
    ap.add_argument("--out-csv", required=True, help="Output summed counts CSV.")
    ap.add_argument("--svg", default=None, help="Optional output SVG path for periodic-table heatmap.")
    ap.add_argument("--title", default="Element Presence per XYZ (file-level)", help="Title for the SVG heatmap.")
    return ap.parse_args()

def load_counts(paths):
    agg = Counter()
    for p in paths:
        with open(p, 'r', encoding='utf-8') as fh:
            r = csv.DictReader(fh)
            if "element" not in r.fieldnames or "count_files_with_element" not in r.fieldnames:
                raise ValueError(f"{p}: missing required headers")
            for row in r:
                el = row["element"].strip()
                try:
                    cnt = int(row["count_files_with_element"])
                except ValueError:
                    try:
                        cnt = float(row["count_files_with_element"])
                    except Exception:
                        continue
                agg[el] += cnt
    return agg

def write_csv(counts, out_csv):
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(["element","count_files_with_element"])
        # include any others (unexpected tokens) at end
        for el in sorted(set(counts.keys()) - set(Z_ORDER)):
            w.writerow([el, counts[el]])
        # include any others (unexpected tokens) at end
        for el in sorted(set(counts.keys()) - set(Z_ORDER)):
            w.writerow([el, counts[el]])

def lerp(a,b,t): return a + (b-a)*t

def rgb(r,g,b): return f'rgb({r},{g},{b})'

def color_for(v, vmax):
    # sequential light→dark blue; keep very light base so labels stay readable
    t = 0.0 if vmax <= 0 else math.log1p(v) / math.log1p(vmax)  # smooth-ish
    r = int(lerp(245,  10, t))
    g = int(lerp(250,  70, t))
    b = int(lerp(255, 160, t))
    return (r, g, b)

def rel_luminance(r,g,b):
    # sRGB → linear → luminance
    def to_lin(c):
        c = c/255.0
        return c/12.92 if c <= 0.04045 else ((c+0.055)/1.055)**2.4
    R, G, B = map(to_lin, (r,g,b))
    return 0.2126*R + 0.7152*G + 0.0722*B

def pick_text_color(r,g,b):
    # white text on dark cells, near-black on light cells
    return "#ffffff" if rel_luminance(r,g,b) < 0.5 else "#111111"

def write_svg_heatmap(counts, out_svg, title):
    # Layout params
    cell = 40
    pad = 12
    label_fs = 12
    title_fs = 16

    # Determine color scale from counts (log-ish for dynamic range)
    vals = [v for v in counts.values() if v > 0]
    vmax = max(vals) if vals else 1
    # log scale mapping t in [0,1]
    def to_t(v):
        if v <= 0: return 0.0
        return math.log10(1+9*v/max(1,vmax))  # compress high values a bit

    # Compute canvas size (18 groups wide, 9 rows to leave ln/an)
    cols = 18
    rows = 9
    width = pad*2 + cols*cell
    height = pad*3 + rows*cell + 40  # extra for title

    def rect(x,y,w,h,fill,stroke="#222"):
        return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="{stroke}" stroke-width="1"/>'

    def text(x,y,txt,anchor="middle",fs=label_fs, weight="500", fill="#111"):
        return (f'<text x="{x}" y="{y}" font-family="Inter,Segoe UI,Arial" '
                f'font-size="{fs}" font-weight="{weight}" text-anchor="{anchor}" '
                f'fill="{fill}">{txt}</text>')
    #def text(x,y,txt,anchor="middle",fs=label_fs, weight="500"):
    #    return f'<text x="{x}" y="{y}" font-family="Inter,Segoe UI,Arial" font-size="{fs}" font-weight="{weight}" text-anchor="{anchor}">{txt}</text>'

    # Build cells
    parts = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    parts.append(text(width/2, pad+title_fs, title, fs=title_fs, weight="600"))

    # Draw grid and labels for all symbols we have positions for
    y0 = pad*2 + title_fs
    for el, (grp, per) in PT_POS.items():
        x = pad + (grp-1)*cell
        y = y0 + (per-1)*cell
        v = counts.get(el, 0)
        r, g, b = color_for(v, vmax)
        fill = rgb(r, g, b)
        fg = pick_text_color(r, g, b)

        parts.append(rect(x, y, cell, cell, fill))
        parts.append(text(x+cell*0.16, y+cell*0.32, el, anchor="start", fs=label_fs, fill=fg))
        if v > 0:
            parts.append(text(x+cell*0.5, y+cell*0.78, str(v), anchor="middle", fs=label_fs-2, weight="400", fill=fg))
    # Simple legend (three ticks)
    lx, ly = pad, y0 + rows*cell + 20
    for i, frac in enumerate([0.0, 0.5, 1.0]):
        v = int(round(vmax*frac))
        cx = lx + i*60
        r, g, b = color_for(v, vmax)      # pass vmax and get RGB tuple
        parts.append(rect(cx, ly, 30, 12, rgb(r,g,b)))
        #parts.append(text(cx+15, ly+28, str(v), fs=label_fs-2))
        parts.append(text(cx+15, ly+28, str(v), fs=label_fs-2, fill="#111"))

    parts.append("</svg>")

    os.makedirs(os.path.dirname(out_svg) or ".", exist_ok=True)
    with open(out_svg, 'w', encoding='utf-8') as fh:
        fh.write("\n".join(parts))

def main():
    args = parse_args()
    csvs = sorted(glob.glob(args.csv_glob))
    if not csvs:
        raise SystemExit(f"No files matched: {args.csv_glob}")

    counts = load_counts(csvs)
    write_csv(counts, args.out_csv)
    print(f"[ok] wrote {args.out_csv}")

    if args.svg:
        write_svg_heatmap(counts, args.svg, args.title)
        print(f"[ok] wrote {args.svg}")

if __name__ == "__main__":
    main()
