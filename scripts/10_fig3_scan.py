"""Figure 3: the saturation-scan heatmap (4 positions x 19 substitutions,
raw predicted ddG relative to GPRP) plus the size-corrected residual at P3,
the one position where the exploratory residual analysis (Results) shows a
pattern distinct from the dominant size trend. Both panels are explicitly
captioned as unvalidated / exploratory, matching the text."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
                     "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "figure.facecolor": "white"})

AA_ORDER = list("GAVLIPFWMSTCYNQDEKRH")  # small->large-ish, grouped by chemistry; WT excluded per-row
d = pd.read_csv(os.path.join(ROOT, "data/results/scan_ddG.csv"))
WT = "GPRP"

fig = plt.figure(figsize=(7.6, 4.6), dpi=300)
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.55)
ax = fig.add_subplot(gs[0])

mat = np.full((4, len(AA_ORDER)), np.nan)
for _, row in d.iterrows():
    j = AA_ORDER.index(row.aa)
    mat[row.position - 1, j] = row.ddG
vmax = np.nanmax(np.abs(mat))
im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(AA_ORDER)), AA_ORDER, fontsize=7.5)
ax.set_yticks(range(4), [f"P{i+1}\n({WT[i]} in GPRP)" for i in range(4)], fontsize=7)
for i in range(4):
    for j, aa in enumerate(AA_ORDER):
        if aa == WT[i]:
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=INK, lw=1.4))
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
cb.set_label(r"predicted $\Delta\Delta G$ vs GPRP (kcal/mol)", fontsize=7.5)
cb.ax.tick_params(labelsize=6.5)
ax.set_title("A  Saturation scan, raw $\\Delta\\Delta G$ -- confounded by side-chain size/charge "
             "(Fig. 2B); not a validated ranking", fontsize=8.3, loc="left", fontweight="bold")
ax.set_xlabel("substituted residue (one-letter code)", fontsize=7.5)

ax2 = fig.add_subplot(gs[1])
r = pd.read_csv(os.path.join(ROOT, "data/results/scan_ddG_residuals.csv"))
r = r.sort_values("resid")
colors = ["#eb6834" if p == 3 else "#2a78d6" for p in r.position]
ax2.axhline(0, color=INK2, lw=0.8)
ax2.bar(range(len(r)), r.resid, color=colors, width=0.82)
ax2.set_xticks(range(len(r)), r.sequence, rotation=90, fontsize=4.6)
ax2.set_ylabel("residual after removing the\nsize trend (kcal/mol)", fontsize=7.3)
ax2.set_title("B  Exploratory only (not pre-registered): $\\Delta\\Delta G$ residual after regressing out "
             "side-chain size (Fig. 2B) -- P3 (Arg→X) substitutions cluster worst",
             fontsize=8.3, loc="left", fontweight="bold")
from matplotlib.patches import Patch
ax2.legend(handles=[Patch(color="#eb6834", label="P3 (native Arg)"),
                    Patch(color="#2a78d6", label="P1/P2/P4")],
          loc="upper left", fontsize=6.5, frameon=False, ncol=2)
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)

fig.savefig(os.path.join(ROOT, "figures/fig3_scan.png"), bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(ROOT, "figures/fig3_scan.pdf"), bbox_inches="tight", facecolor="white")
print("saved fig3")
