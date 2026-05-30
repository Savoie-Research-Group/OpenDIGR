#!/usr/bin/env python3
"""
element_heatmap_grouped.py

Periodic-table heatmap where each cell is colored by element *group*
(same palette as element_group_proportion_stripe.py) and the color
intensity within a group scales with `count_files_with_element` from
hist_input.csv (log-scaled).

Output: panel_a_heatmap_v2.{svg,png}
"""
import csv
import math
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgb

HERE = os.path.dirname(os.path.abspath(__file__))
IN_CSV = os.path.join(HERE, "hist_input.csv")
OUT_SVG = os.path.join(HERE, "panel_a_heatmap_v2.svg")
OUT_PNG = os.path.join(HERE, "panel_a_heatmap_v2.png")

# ---------------------------------------------------------------------------
# Group definitions and palette — keep in lock-step with proportion_stripe.py
# ---------------------------------------------------------------------------

GROUP_COLOR = {
    "ALKALI":           "#F48FB1",                          # pink
    "ALKALINE_EARTH":   "#F6D365",                          # gold
    "TRANSITION":       "#66BB6A",                          # green
    "POST_TRANSITION":  "#4DD0E1",                          # cyan
    "METALLOIDS":       "#84CDEE",                          # light blue (132,205,238)
    "LANTHANIDES":      "#8a2be2",                          # blue-violet (incl. La)
    # Extra categories present in hist_input but not in the stripe:
    "NONMETAL":         "#84CDEE",                          # same shade as metalloids (Figure-2 convention)
    "NOBLE_GAS":        "#E0E0E0",                          # very light gray
    "HYDROGEN":         "#E0E0E0",                          # light gray (often shown apart)
}

GROUPS = {
    "ALKALI":          ["Li", "Na", "K", "Rb", "Cs", "Fr"],
    "ALKALINE_EARTH":  ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"],
    "TRANSITION":      ["Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                        "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
                        "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg"],
    "POST_TRANSITION": ["Al","Ga","In","Sn","Tl","Pb","Bi","Po"],
    "METALLOIDS":      ["B","Si","Ge","As","Sb","Te"],
    "LANTHANIDES":     ["La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
                        "Ho","Er","Tm","Yb","Lu"],
    "NONMETAL":        ["C","N","O","F","P","S","Cl","Se","Br","I","At"],
    "NOBLE_GAS":       ["He","Ne","Ar","Kr","Xe","Rn"],
    "HYDROGEN":        ["H"],
}

EL_GROUP = {el: g for g, els in GROUPS.items() for el in els}

# ---------------------------------------------------------------------------
# Periodic-table layout: (row, col) in a standard 7-row + 1-row-Ln grid.
# Cols 1..18.
# ---------------------------------------------------------------------------

LAYOUT = {
    "H":  (1, 1),  "He": (1, 18),
    "Li": (2, 1),  "Be": (2, 2),
    "B":  (2, 13), "C": (2, 14), "N": (2, 15), "O": (2, 16), "F": (2, 17), "Ne": (2, 18),
    "Na": (3, 1),  "Mg": (3, 2),
    "Al": (3, 13), "Si": (3, 14), "P": (3, 15), "S": (3, 16), "Cl": (3, 17), "Ar": (3, 18),
    "K":  (4, 1),  "Ca": (4, 2),
    "Sc": (4, 3),  "Ti": (4, 4), "V": (4, 5), "Cr": (4, 6), "Mn": (4, 7),
    "Fe": (4, 8),  "Co": (4, 9), "Ni": (4, 10), "Cu": (4, 11), "Zn": (4, 12),
    "Ga": (4, 13), "Ge": (4, 14), "As": (4, 15), "Se": (4, 16), "Br": (4, 17), "Kr": (4, 18),
    "Rb": (5, 1),  "Sr": (5, 2),
    "Y":  (5, 3),  "Zr": (5, 4), "Nb": (5, 5), "Mo": (5, 6), "Tc": (5, 7),
    "Ru": (5, 8),  "Rh": (5, 9), "Pd": (5, 10), "Ag": (5, 11), "Cd": (5, 12),
    "In": (5, 13), "Sn": (5, 14), "Sb": (5, 15), "Te": (5, 16), "I": (5, 17), "Xe": (5, 18),
    "Cs": (6, 1),  "Ba": (6, 2),  "La": (6, 3),
    "Hf": (6, 4),  "Ta": (6, 5),  "W": (6, 6), "Re": (6, 7),
    "Os": (6, 8),  "Ir": (6, 9), "Pt": (6, 10), "Au": (6, 11), "Hg": (6, 12),
    "Tl": (6, 13), "Pb": (6, 14), "Bi": (6, 15), "Po": (6, 16), "At": (6, 17), "Rn": (6, 18),
    # Lanthanide row, displayed below main table (row 8 here, cols 4..17)
    "Ce": (8, 4),  "Pr": (8, 5),  "Nd": (8, 6), "Pm": (8, 7), "Sm": (8, 8),
    "Eu": (8, 9),  "Gd": (8, 10), "Tb": (8, 11), "Dy": (8, 12),
    "Ho": (8, 13), "Er": (8, 14), "Tm": (8, 15), "Yb": (8, 16), "Lu": (8, 17),
}

# ---------------------------------------------------------------------------
# Load counts
# ---------------------------------------------------------------------------

counts = {}
with open(IN_CSV) as fh:
    for row in csv.DictReader(fh):
        counts[row["element"]] = int(row["count_files_with_element"])

# Log scale for intensity (counts span 2 -> 330k)
log_min = math.log(2)
log_max = max(math.log(c + 1) for c in counts.values())

def intensity(el):
    c = counts.get(el)
    if c is None or c < 1:
        return 0.0
    return (math.log(c + 1) - log_min) / (log_max - log_min)

def mix_white(rgb_tuple, t):
    """Interpolate between white (t=0) and the saturated color (t=1)."""
    t = max(0.05, min(1.0, t))  # floor so no element is pure white
    r, g, b = rgb_tuple
    return (1 - t + t * r,
            1 - t + t * g,
            1 - t + t * b)

# ---------------------------------------------------------------------------
# Draw
# ---------------------------------------------------------------------------

CELL = 1.7                  # only slightly bigger than v3; text will dominate the cell
GAP  = 0.05
N_COLS = 18
N_ROWS = 8                  # rows 1..8 (8 = lanthanide row)
LEGEND_H = 0.2              # small bottom margin; the color legend has been removed

fig_w = N_COLS * CELL + 1.2
fig_h = N_ROWS * CELL + LEGEND_H + 0.6
fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=180)
ax.set_xlim(0, N_COLS * CELL + 0.4)
ax.set_ylim(0, N_ROWS * CELL + LEGEND_H + 0.4)
ax.set_aspect("equal")
ax.axis("off")

# Draw each cell that has a layout position. All cell y values are shifted
# up by LEGEND_H so the bottom strip (y < LEGEND_H) is reserved for the legend.
for el, (row, col) in LAYOUT.items():
    grp = EL_GROUP.get(el)
    if grp is None:
        continue
    base = to_rgb(GROUP_COLOR[grp])
    t = intensity(el)
    fc = mix_white(base, t)
    x = (col - 1) * CELL + 0.2
    # Invert row so row 1 is at top; add LEGEND_H to lift above legend
    y = (N_ROWS - row) * CELL + LEGEND_H + 0.2
    rect = Rectangle((x + GAP/2, y + GAP/2), CELL - GAP, CELL - GAP,
                     facecolor=fc, edgecolor="#222", linewidth=0.4)
    ax.add_patch(rect)
    # Element symbol
    ax.text(x + CELL/2, y + CELL - 0.34, el,
            ha="center", va="center", fontsize=32, fontweight="bold",
            color="#111")
    # Count (explicit, no comma, bold). Sized so the 6-digit max
    # (330446 for H) still fits within the cell.
    c = counts.get(el)
    if c is not None:
        ax.text(x + CELL/2, y + 0.46, str(c),
                ha="center", va="center", fontsize=23,
                fontweight="bold", color="#111")

plt.tight_layout()
plt.savefig(OUT_SVG, bbox_inches="tight")
plt.savefig(OUT_PNG, bbox_inches="tight")
print(f"Wrote: {OUT_SVG}", file=sys.stderr)
print(f"Wrote: {OUT_PNG}", file=sys.stderr)
