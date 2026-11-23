#!/usr/bin/env python
"""Figure S11 — Annotation-independent classification separates folial Bergmann-like ML astroglia
from boundary-associated Astro-NT cells (review 3e/3f, design Antonin, 2026-08-17).

(a) design schematic — training on anatomically unambiguous references only, 0 ML cells;
(b) BG-score vs Astro-NT-score: references + the 519 ML astroglia projected with labels withheld;
(c) concordance of the eight classifications (% BG / % Astro-NT per census group);
(d) majority predictions mapped on section ABCA-3.005 — predicted Astro-NT concentrate at borders.
Reads: annotfree_scores_cells.csv, annotfree_summary.csv (build_annotation_free_identity.py).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib import patheffects as pe

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, ENTITY, panel_letter, save, W2
from geom_layers import layer_regions, paint_layers

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"
ZH = BASE + "/datasets/zhang_zhuang_2023_merfish"
MAP_SEC = "Zhuang-ABCA-3.005"
CB_COL = ENTITY["Astro_CB"]; NT_COL = ENTITY["Astro_NT"]; BG_COL = GROUP["BG"]
LAYER_TINT = {"ML": "#dcebf9", "PCL": "#f6d9a8", "GCL": "#e7ddf3", "WM": "#e4e7eb"}

cells = pd.read_csv(TAB + "/annotfree_scores_cells.csv")
summary = pd.read_csv(TAB + "/annotfree_summary.csv")
ml = cells[cells["set"] == "ML-census"]
assert len(ml) == 519, len(ml)
assert (cells["set"] == "ref BG-PCL").sum() == 800

fig = plt.figure(figsize=(W2, 7.4))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1.0, 1.15], height_ratios=[1.0, 0.95],
              wspace=0.34, hspace=0.36, top=0.95, bottom=0.06, left=0.075, right=0.985)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])

# ---------- a: design schematic ----------
axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.axis("off")
def box(x, y, w, h, text, fc, fs=7.0, tc="#111111"):
    axA.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                                 facecolor=fc, edgecolor="#555555", lw=0.7, clip_on=False))
    axA.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc)
# box widths sized to their text (review O.G. 2026-10-23: labels overflowed the boxes)
box(0.03, 0.80, 0.50, 0.13, "BG reference\nsomata in the PCL (n=800)", "#d9e6f7")
box(0.03, 0.62, 0.50, 0.13, "Astro-NT reference\ninside the WM (n=800)", "#fbeccc")
box(0.03, 0.44, 0.50, 0.13, "Astro-CB reference\ninside the GCL (n=800)", "#d9f0e6")
box(0.58, 0.57, 0.39, 0.22, "4 classifiers\n× 2 gene spaces\n(full / contam.-free)", "#efefef")
box(0.20, 0.12, 0.64, 0.16, "519 ML astroglia\nlabels withheld → classified", "#f4d7d7")
for y0 in (0.865, 0.685, 0.505):
    axA.add_patch(FancyArrowPatch((0.535, y0), (0.575, 0.70), arrowstyle="-|>",
                                  mutation_scale=9, color="#555555", lw=0.9, clip_on=False))
axA.add_patch(FancyArrowPatch((0.78, 0.565), (0.62, 0.29), arrowstyle="-|>",
                              mutation_scale=10, color="#555555", lw=1.0, clip_on=False))
axA.text(0.52, 0.045, "0 ML cells used for training", fontsize=7.8, fontweight="bold",
         color="#8a2525", ha="center", va="center")
panel_letter(axA, "a", dx=0.0, dy=1.02)

# ---------- b: scatter, ML cells projected without labels ----------
for name, col in [("ref BG-PCL", BG_COL), ("ref AstroNT-WM", NT_COL), ("ref CB-GCL", CB_COL)]:
    s = cells[cells["set"] == name]
    axB.scatter(s.bg_score, s.nt_score, s=7, c=col, alpha=0.22, linewidths=0)
axB.scatter(ml.bg_score, ml.nt_score, s=13, facecolors="none", edgecolors="#111111",
            linewidths=0.6, alpha=0.85, zorder=6)
axB.axhline(0, c="#cccccc", lw=0.75, zorder=0); axB.axvline(0, c="#cccccc", lw=0.75, zorder=0)
axB.set_xlabel("BG signature score\n(BG-PCL vs Astro-NT-WM; no ML cell in either signature)")
axB.set_ylabel("Astro-NT signature score")
h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=BG_COL, markersize=6, label="BG-PCL reference"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=NT_COL, markersize=6, label="Astro-NT-WM reference"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=CB_COL, markersize=6, label="Astro-CB-GCL reference"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#111111",
            markersize=6, label="ML astroglia (labels withheld)")]
axB.legend(handles=h, fontsize=7.0, loc="upper right", frameon=False, handletextpad=0.3)
panel_letter(axB, "b", dx=-0.20, dy=1.02)

# ---------- c: concordance of the eight classifications ----------
GRPS = ["ML-Bergmann", "ML-NT-border", "ML-NT-interior", "ML-CB"]
RUNS = [f"{m}_{t}" for t in ("full", "clean") for m in ("logreg", "rforest", "kNN15", "corrcent")]
def mat(assigned):
    M = np.zeros((len(GRPS), len(RUNS)))
    for i, g in enumerate(GRPS):
        for k, run in enumerate(RUNS):
            v = summary[(summary.group == g) & (summary.run == run) & (summary.assigned == assigned)].pct
            M[i, k] = float(v.iloc[0]) if len(v) else np.nan
    return M
axC.axis("off")
ns = {g: int(summary[summary.group == g].n.iloc[0]) for g in GRPS}
for bi, (assigned, cmap, y0) in enumerate([("BG", "Blues", 0.56), ("Astro-NT", "Oranges", 0.02)]):
    M = mat(assigned)
    axc = axC.inset_axes([0.14, y0, 0.86, 0.34])
    axc.imshow(M, cmap=cmap, vmin=0, vmax=100, aspect="auto")
    for i in range(M.shape[0]):
        for k in range(M.shape[1]):
            axc.text(k, i, f"{M[i, k]:.0f}", ha="center", va="center", fontsize=7.0,
                     color="white" if M[i, k] > 60 else "#222222")
    axc.set_yticks(range(len(GRPS)))
    axc.set_yticklabels([f"{g} ({ns[g]})" for g in GRPS], fontsize=7.0)
    if bi == 1:
        axc.set_xticks(range(len(RUNS)))
        axc.set_xticklabels([r.replace("_", "\n") for r in RUNS], fontsize=7.0,
                            rotation=45, ha="right", rotation_mode="anchor")
    else:
        axc.set_xticks([])
    axc.set_title(f"% classified as {assigned}", fontsize=7.6, loc="left", pad=2)
    axc.tick_params(length=0)
panel_letter(axC, "c", dx=0.0, dy=1.02)

# ---------- d: majority predictions on the map ----------
mv = pd.read_parquet(ZH + "/cb_cells_multivote.parquet")
zz = mv[mv.brain_section_label == MAP_SEC]
reg = layer_regions(zz.x, zz.y, zz.is_gran, zz.subclass.astype(str).str.contains("Oligo"), umpp=1000.0)
paint_layers(axD, reg, colors=LAYER_TINT, pcl_line="#b8560f", pcl_lw=0.8)
sm = ml[ml.section == MAP_SEC]
PRED_COL = {"BG": BG_COL, "Astro-NT": NT_COL, "Astro-CB": CB_COL}
for pcls, col in PRED_COL.items():
    s = sm[sm.pred_majority == pcls]
    if pcls == "BG":
        axD.scatter(s.x, s.y, s=46, marker="*", c=col, edgecolors="#1c3a5e", linewidths=0.4, zorder=6)
    else:
        axD.scatter(s.x, s.y, s=20, c=col, edgecolors="#333333", linewidths=0.4, zorder=6)
axD.set_aspect("equal"); axD.axis("off")
axD.set_title("majority prediction, section ABCA-3.005", fontsize=8, loc="right", pad=2)
h2 = [Line2D([0], [0], marker="*", color="w", markerfacecolor=BG_COL, markeredgecolor="#1c3a5e",
             markersize=10, label="predicted BG"),
      Line2D([0], [0], marker="o", color="w", markerfacecolor=NT_COL, markeredgecolor="#333333",
             markersize=6, label="predicted Astro-NT"),
      Line2D([0], [0], marker="o", color="w", markerfacecolor=CB_COL, markeredgecolor="#333333",
             markersize=6, label="predicted Astro-CB")]
axD.legend(handles=h2, fontsize=7.0, loc="lower left", frameon=False, handletextpad=0.3)
panel_letter(axD, "d", dx=0.0, dy=1.02)

save(fig, "FigureS11_annotation_free")
