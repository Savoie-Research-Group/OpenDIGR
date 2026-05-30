#
# proportion_stripe.py
# Create a thin, long composition stripe for element-group shares.

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe

# ---- Data (replace here if needed) ----
counts = {
    "TRANSITION":       179401,
    "METALLOIDS":        28454,
    "ALKALI":            10370,
    "POST_TRANSITION":    7642,
    "ALKALINE_EARTH":     4200,
    "LANTHANIDES":         584,
}

# ---- Color scheme (match your heatmap) ----
# Alkali: pink, Alkaline earth: yellow, Transition: green,
# Post-transition: cyan, Metalloid: light blue (slightly bluer),
# Lanthanides: grey
colors = {
    "ALKALI":           "#F48FB1",                      # pink
    "ALKALINE_EARTH":   "#F6D365",                      # yellow
    "TRANSITION":       "#66BB6A",                      # green
    "POST_TRANSITION":  "#4DD0E1",                      # cyan
    "METALLOIDS":       (132/255, 205/255, 238/255),    # light blue (a bit bluer)
    "LANTHANIDES":      "#BFBFBF",                      # grey
}

# ---- Compute percentages & order ----
total = sum(counts.values())
items = sorted(((k, v, v/total) for k, v in counts.items()),
               key=lambda x: x[2], reverse=True)

# ---- Figure: thin & long ----
fig_w, fig_h = 10, 0.6  # inches; tweak width to match your heatmap width
fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# ---- Draw stripe ----
x = 0.0
pad = 0.001  # tiny gap between segments; set 0 for no gaps
for name, count, frac in items:
    w = frac - pad if frac > pad else frac
    rect = Rectangle((x, 0), w, 1, facecolor=colors[name], edgecolor="none")
    ax.add_patch(rect)

    # Direct label (bold, black). If very small, print outside above the bar.
    label = f"{name.replace('_',' ')} {100*frac:.1f}%"
    text_effects = [pe.withStroke(linewidth=1.2, foreground="white")]
    if w >= 0.06:  # label inside if there is room
        ax.text(x + w/2, 0.5, label,
                ha="center", va="center", fontsize=10, color="black",
                fontweight="bold", path_effects=text_effects)
    else:
        ax.text(min(x + w + 0.003, 0.995), 0.5, label,
                ha="left", va="center", fontsize=9, color="black",
                fontweight="bold")

    x += frac

# ---- Save ----
plt.savefig("element_group_proportion_stripe.svg", bbox_inches="tight")
plt.savefig("element_group_proportion_stripe.png", bbox_inches="tight")
plt.close(fig)

print("Saved: element_group_proportion_stripe.svg and .png")
