"""Figure 1: (A) rendered hole-a pocket with GPRP (PyMOL image), (B) contact map of the
template peptide residues vs receptor residues within 4.2 A (distances from 2FIB)."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import knobhole_lib as kl  # noqa: E402

ROOT = os.path.dirname(HERE)
INK, INK2 = "#0b0b0b", "#52514e"
ACIDIC, POLAR, AROM = "#c0524f", "#2a8a8a", "#a07a2a"      # text colours echo the 3-D colouring
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2,
                     "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})

rec, ca, pep = kl.read_template()
names = ["Gly1", "Pro2", "Arg3", "Pro4"]
pep_atoms = [np.array([x for n, (x, e) in a.items() if e != "H"]) for _, _, a in pep]
cols = []
for seq, (resn, atoms) in sorted(rec.items()):
    xyz = kl._heavy_xyz(atoms)
    d = [np.linalg.norm(xyz[:, None] - p[None], axis=2).min() for p in pep_atoms]
    if min(d) <= 4.2:
        cols.append((seq, resn.title(), d))
M = np.array([c[2] for c in cols]).T
labels = [f"{c[1]}{c[0]}" for c in cols]

fig = plt.figure(figsize=(7.2, 3.9), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.08)
axa = fig.add_subplot(gs[0])
axa.imshow(plt.imread(os.path.join(ROOT, "figures/src/fig1a_pocket.png")))
axa.axis("off")
axa.text(0.0, 1.02, "A", transform=axa.transAxes, fontsize=12, fontweight="bold", va="bottom")
axa.text(0.5, 1.02, "Hole a with bound GPRP (orange), PDB 2FIB", transform=axa.transAxes,
         ha="center", va="bottom", fontsize=8.5, color=INK2)
axa.text(0.03, 0.04, "green sphere: structural Ca$^{2+}$\nsalmon: Asp   teal: Gln/Lys/His/Arg\nsand: Phe/Thr",
         transform=axa.transAxes, fontsize=6.8, color=INK2, va="bottom")

axb = fig.add_subplot(gs[1])
cmap = LinearSegmentedColormap.from_list("blu", ["#0d366b", "#2a78d6", "#9ec5f4", "#f7f7f5"])
im = axb.imshow(np.clip(M, 2.5, 4.2), cmap=cmap, aspect="auto", vmin=2.5, vmax=4.2)
axb.set_yticks(range(4), names)
axb.set_xticks(range(len(labels)), labels, rotation=90, fontsize=7)
for tick, lab in zip(axb.get_xticklabels(), labels):
    num = int(''.join(ch for ch in lab if ch.isdigit()))
    tick.set_color(ACIDIC if num in (364, 330, 297, 301) else POLAR if num in (329, 338, 340, 375, 339) else
                   AROM if num in (295, 322, 305) else INK2)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if M[i, j] <= 4.2:
            axb.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=5.8,
                     color="white" if M[i, j] < 3.3 else INK)
axb.set_xticks(np.arange(-.5, len(labels), 1), minor=True)
axb.set_yticks(np.arange(-.5, 4, 1), minor=True)
axb.grid(which="minor", color="white", linewidth=1.2)
axb.tick_params(which="minor", length=0)
for s in axb.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=axb, fraction=0.04, pad=0.02)
cb.set_label("closest heavy-atom distance (Å)", fontsize=7.5)
cb.ax.tick_params(labelsize=7)
axb.text(0.0, 1.02, "B", transform=axb.transAxes, fontsize=12, fontweight="bold", va="bottom")
axb.text(0.5, 1.02, "Peptide-receptor contacts ≤ 4.2 Å", transform=axb.transAxes, ha="center",
         va="bottom", fontsize=8.5, color=INK2)
os.makedirs(os.path.join(ROOT, "figures"), exist_ok=True)
fig.savefig(os.path.join(ROOT, "figures/fig1_pocket.png"), bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(ROOT, "figures/fig1_pocket.pdf"), bbox_inches="tight", facecolor="white")
print("saved fig1;", len(labels), "contact residues:", labels)
