#!/usr/bin/env python
"""Figure S10 — sensitivity of the main results to the FC selection margin and to manual review.

Reads tables/selection_sensitivity.csv (build_selection_sensitivity.py). Rows a-h: MERFISH,
four animals; every readout recomputed with the canonical definition for each cohort:
  validated   — manually validated FC at margins 35-100 µm (45 = the cohort of the study; the
                35-µm cohort adds the candidates of the 35-45 µm band accepted in a dedicated review);
  geometric   — every geometric candidate at the margin, no manual review (35-100 µm);
  rejected    — the candidates discarded by manual review (45 µm), on their own.
Row i-k: validated subsets on the two other platforms (Stereo-seq: depth above the PCL;
Xenium: distance to the nearest Purkinje soma, canonical rule >= 30 µm).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, INK, W2, panel_letter, save
from project_paths import TABLES

use_style()
t = pd.read_csv(TABLES / "selection_sensitivity.csv")
m = t[t.dataset == "MERFISH"]
val = m[m.arm == "validated"].sort_values("margin_um")
geo = m[m.arm == "geometric"].sort_values("margin_um")
rej = m[m.arm == "rejected"].iloc[0]
canon = val[val.margin_um == 45].iloc[0]

C_VAL, C_GEO, C_REF = GROUP["FC"], "#4d4d4d", "#8a8a8a"
XT = [35, 45, 55, 65, 80, 100]


def draw(ax, col, ylabel, ylim, ref=None, ref_label=None, sig_col=None):
    """One MERFISH readout against the margin, all cohorts."""
    if ref is not None:
        ax.axhline(ref, color=C_REF, lw=0.9, ls=(0, (4, 3)), zorder=1)
        if ref_label:
            ax.annotate(ref_label, (105, ref), xytext=(0, 2.5), textcoords="offset points",
                        fontsize=7, color=C_REF, ha="right", va="bottom")
    ax.plot(geo.margin_um, geo[col], "-", color=C_GEO, lw=1.2, zorder=2)
    ax.plot(val.margin_um, val[col], "-", color=C_VAL, lw=1.6, zorder=4)
    for df, color, marker, z in [(geo, C_GEO, "s", 2), (val, C_VAL, "o", 4)]:
        for _, r in df.iterrows():
            sig = True if sig_col is None else (r[sig_col] < 0.05)
            ax.plot(r.margin_um, r[col], marker, ms=4.5, mec=color, mew=1.0,
                    mfc=color if sig else "white", zorder=z + 1)
    sig = True if sig_col is None else (rej[sig_col] < 0.05)
    ax.plot(45, rej[col], "X", ms=6, mec=C_GEO, mew=0.9, mfc=C_GEO if sig else "white", zorder=3)
    ax.axvline(45, color=INK if isinstance(INK, str) else "#222222", lw=0.75, ls=":", zorder=0)
    ax.set_xticks(XT)
    ax.set_xticklabels([str(x) for x in XT], fontsize=7)
    ax.set_xlim(29, 106)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)


fig = plt.figure(figsize=(W2, 8.2))
gs = fig.add_gridspec(3, 4, wspace=0.78, hspace=0.85, left=0.085, right=0.985, top=0.875, bottom=0.06)
ax = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(4)]

draw(ax[0], "n_FC", "FC in cohort (n)", (0, 1050))
draw(ax[1], "knn_pct_BG", "kNN: % BG neighbors", (60, 103), ref=canon.knn_base_pct_BG, ref_label="chance")
draw(ax[2], "centroid_ratio", "centroid ratio", (0, 5.8), ref=1, ref_label="equidistant")
draw(ax[3], "silhouette_FC_BG", "silhouette FC | BG", (-0.06, 0.36), ref=canon.silhouette_BG_par,
     ref_label="BG | parench.")
draw(ax[4], "med_abs_d_FC_BG", "median |Δlog2| FC–BG", (0, 0.145), ref=canon.med_abs_d_BG_par,
     ref_label="BG–parench.")
draw(ax[5], "pct_genes_gt025_FC_BG", "% genes |Δlog2| > 0.25", (0, 36), ref=canon.pct_genes_gt025_BG_par,
     ref_label="BG–parench.")
draw(ax[6], "n_DE_hits", "matched DE hits (n)", (-0.6, 12))
ax[6].set_yticks([0, 3, 6, 9, 12])
draw(ax[7], "lean_delta", "astrocytic shift Δ", (-0.06, 0.32), ref=0, sig_col="lean_p")
for a in ax[4:]:
    a.set_xlabel("selection margin (µm)")
for i, a in enumerate(ax):
    panel_letter(a, "abcdefgh"[i], dx=-0.36, dy=1.15)

handles = [
    Line2D([], [], color=C_VAL, marker="o", ms=4.5, lw=1.6, label="manually validated FC (45 µm = cohort of the study)"),
    Line2D([], [], color=C_GEO, marker="s", ms=4.5, lw=1.2, label="all geometric candidates, no manual review"),
    Line2D([], [], color=C_GEO, marker="X", ms=6, lw=0, label="candidates rejected by manual review"),
    Line2D([], [], color=C_REF, lw=0.9, ls=(0, (4, 3)), label="reference level"),
]
fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=7,
           bbox_to_anchor=(0.53, 0.995), columnspacing=1.4, handlelength=2.6)
fig.text(0.5, ax[0].get_position().y1 + 0.038, "MERFISH (four animals)", fontsize=9, fontweight="bold",
         ha="center", va="bottom")

# ---------------- other platforms: validated subsets ----------------
OTH = {"Stereo-seq": ("#7a5195", "o"), "Xenium": ("#ef8a17", "^")}
specs = [("knn_pct_BG", "kNN: % BG neighbors", (40, 103), None, None),
         ("centroid_ratio", "centroid ratio", (0, 2.6), 1, None),
         ("lean_delta", "astrocytic shift Δ", (-0.12, 0.36), 0, "lean_p")]
axo = [fig.add_subplot(gs[2, j]) for j in range(3)]
for a, (col, lab, yl, ref, sig_col) in zip(axo, specs):
    if ref is not None:
        a.axhline(ref, color=C_REF, lw=0.9, ls=(0, (4, 3)), zorder=1)
    for name, (color, marker) in OTH.items():
        d = t[(t.dataset == name) & (t.arm == "validated")].sort_values("margin_um")
        a.plot(d.margin_um, d[col], "-", color=color, lw=1.4, zorder=2)
        for _, r in d.iterrows():
            sig = True if sig_col is None else (r[sig_col] < 0.05)
            a.plot(r.margin_um, r[col], marker, ms=4.8, mec=color, mew=1.0,
                   mfc=color if sig else "white", zorder=3)
        if col == "knn_pct_BG":
            a.plot([d.margin_um.min(), d.margin_um.max()], [d.knn_base_pct_BG.iloc[0]] * 2,
                   color=color, lw=0.9, ls=(0, (4, 3)), zorder=1)
    a.set_xticks([30, 45, 65, 80, 100])
    a.set_xticklabels(["30", "45", "65", "80", "100"], fontsize=7)
    a.set_xlim(24, 106)
    a.set_ylim(*yl)
    a.set_ylabel(lab)
    a.set_xlabel("selection margin (µm)")
for i, a in enumerate(axo):
    panel_letter(a, "ijk"[i], dx=-0.36, dy=1.15)
axl = fig.add_subplot(gs[2, 3])
axl.axis("off")
axl.legend(handles=[Line2D([], [], color=c, marker=mk, ms=4.8, lw=1.4, label=n) for n, (c, mk) in OTH.items()]
           + [Line2D([], [], color=C_REF, marker="o", ms=4.5, lw=0, mfc="white", mec=C_REF,
                     label="open symbol: p ≥ 0.05\n(shift panels h, k)"),
              Line2D([], [], color=C_REF, lw=0.9, ls=(0, (4, 3)), label="chance / null level")],
           loc="center left", frameon=False, fontsize=7, bbox_to_anchor=(-0.25, 0.5), handlelength=2.4)
fig.text(0.5, axo[0].get_position().y1 + 0.038, "Stereo-seq and Xenium (validated FC)", fontsize=9,
         fontweight="bold", ha="center", va="bottom")

save(fig, "FigureS10_threshold")
