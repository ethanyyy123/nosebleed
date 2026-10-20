"""Figure 2: the H1 validation panel. Two axes: (A) bar chart of ddG (mean +/-
SEM) for every H1 comparison, coloured by whether it meets its own
pre-registered criterion; (B) the saturation scan's own pre-registered
size-bias diagnostic (ddG vs change in side-chain heavy-atom count, all 76
scanned variants), since that diagnostic is what ultimately determines
whether the scan in Figure 3 can be read as more than a methods note."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
GOOD, CRIT = "#0ca30c", "#d03b3b"
BLUE = "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
                     "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.facecolor": SURF, "figure.facecolor": "white"})

h1 = pd.read_csv(os.path.join(ROOT, "data/results/h1_panel.csv"))
scan = pd.read_csv(os.path.join(ROOT, "data/results/scan_ddG.csv"))

fig, (axa, axb) = plt.subplots(1, 2, figsize=(8.0, 4.0), dpi=300, gridspec_kw=dict(width_ratios=[1.0, 1.0], wspace=0.40))

# --- A: H1 bars ---
labels = {"GPR vs GPRP (6 A)": "GPR\n(H1a)",
          "D364A vs WT (7 A)": "D364A\n(H1b)", "D364H vs WT (7 A)": "D364H\n(H1b)",
          "D364V vs WT (7 A)": "D364V\n(H1b)", "F293A vs WT (7 A)": "F293A\n(H1c)",
          "S300A vs WT (7 A)": "S300A\n(H1c)", "M336A vs WT (7 A)": "M336A\n(H1c)",
          "N337A vs WT (7 A)": "N337A\n(H1c)"}
h1 = h1.set_index("comparison").loc[list(labels)].reset_index()
x = np.arange(len(h1))
colors = [GOOD if p else CRIT for p in h1.passed]
axa.bar(x, h1.ddG, color=colors, width=0.62, zorder=3)
axa.axhline(0, color=INK2, lw=0.8)
axa.axhspan(-3, 3, color=INK2, alpha=0.08, zorder=0)
axa.set_ylim(-9, 14.5)
for xi, row in zip(x, h1.itertuples()):
    axa.text(xi, row.ddG + (0.6 if row.ddG >= 0 else -0.6), f"n={row.n_test}",
             ha="center", va="bottom" if row.ddG >= 0 else "top", fontsize=6.3, color=INK2)
axa.set_xticks(x, [labels[c] for c in h1.comparison], fontsize=6.6)
axa.set_ylabel(r"$\Delta\Delta G$ vs same-shell WT (kcal/mol)", fontsize=8)
axa.set_title("A  H1 validation panel", fontsize=9, loc="left", fontweight="bold")
axa.text(0.02, 0.97, "shaded band: |$\\Delta\\Delta G$| < 3 kcal/mol (pre-set no-effect threshold)",
         transform=axa.transAxes, ha="left", va="top", fontsize=6.0, color=INK2)
from matplotlib.patches import Patch
axa.legend(handles=[Patch(color=GOOD, label="meets its criterion"),
                    Patch(color=CRIT, label="fails its criterion")],
          loc="lower left", fontsize=6.5, frameon=False)
for s in ("top", "right"):
    axa.spines[s].set_visible(False)

# --- B: size-bias diagnostic ---
d = scan.copy()
d["dsize"] = d.sc_heavy - d.sc_heavy_wt
slope, intercept, r, p, se = stats.linregress(d.dsize, d.ddG)
xs = np.linspace(d.dsize.min() - 0.5, d.dsize.max() + 0.5, 50)
axb.axhline(0, color=INK2, lw=0.8, zorder=1)
axb.scatter(d.dsize, d.ddG, s=16, c=BLUE, alpha=0.75, linewidths=0, zorder=3)
axb.plot(xs, intercept + slope * xs, color=CRIT, lw=1.6, zorder=2)
axb.text(0.03, 0.05, f"$\\rho$ = -0.79, p = 3$\\times10^{{-17}}$ (pooled, n=76)\nslope {slope:.2f} kcal/mol per added heavy atom",
         transform=axb.transAxes, fontsize=6.8, color=INK2, va="bottom")
axb.set_xlabel("change in side-chain heavy-atom count vs GPRP", fontsize=8)
axb.set_ylabel(r"predicted $\Delta\Delta G$ (kcal/mol)", fontsize=8)
axb.set_title("B  pre-registered size-bias diagnostic\n(all 76 scanned variants)", fontsize=9, loc="left", fontweight="bold")
for s in ("top", "right"):
    axb.spines[s].set_visible(False)

fig.savefig(os.path.join(ROOT, "figures/fig2_h1_and_sizebias.png"), bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(ROOT, "figures/fig2_h1_and_sizebias.pdf"), bbox_inches="tight", facecolor="white")
print("saved fig2")
