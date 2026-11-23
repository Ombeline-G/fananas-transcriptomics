#!/usr/bin/env python
"""Supplementary figures built from canonical data:
S6  ONE cross-platform heatmap of all functional + metabolic modules (review\n    point 27; replaces the former S6/S7 pair; complete table:\n    TableS_functional_metabolic_modules.csv).
S8  per-section paired lean (harmonised protocol) — pseudoreplication control.
S10 is no longer built here: see FigureS10_threshold.py (name of this file kept for the README).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, DSL, DS_ORDER, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"
C_NS = "#b8b6b0"

# ================= S6 — one cross-platform module heatmap (review point 27) =================
t = pd.read_csv(TAB + "/TableS_functional_metabolic_modules.csv")
DS_LABELS = ["MERFISH", "Stereo-seq", "Xenium"]
mods = (t.drop_duplicates(["module_class", "module"])
          .sort_values(["module_class", "module"], ascending=[True, True]))
mods = pd.concat([mods[mods.module_class == "functional"],
                  mods[mods.module_class == "metabolic"]])
rows_order = list(mods.module)
piv_d = t.pivot_table(index="module", columns="dataset", values="delta").reindex(rows_order)
piv_q = t.pivot_table(index="module", columns="dataset", values="padj").reindex(rows_order)
piv_n = t.pivot_table(index="module", columns="dataset", values="n_genes").reindex(rows_order)
n_func = int((mods.module_class == "functional").sum())

fig, ax = plt.subplots(figsize=(W2 * 0.62, 5.4))
fig.subplots_adjust(left=0.40, top=0.90, bottom=0.14, right=0.97)
vmax = 0.5
M = piv_d[DS_LABELS].to_numpy(dtype=float)
im = ax.imshow(np.clip(M, -vmax, vmax), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        q = piv_q[DS_LABELS].to_numpy()[i, j]
        n = piv_n[DS_LABELS].to_numpy()[i, j]
        if np.isnan(M[i, j]):
            ax.text(j, i, "n.m.", ha="center", va="center", fontsize=7, color="#898781")
        else:
            mark = "*" if (not np.isnan(q) and q < 0.05) else ""
            small = f" n={int(n)}" if n <= 2 else ""
            ax.text(j, i, f"{M[i, j]:+.2f}{mark}{small}", ha="center", va="center",
                    fontsize=7, color="#0b0b0b" if abs(M[i, j]) < 0.30 else "white")
ax.axhline(n_func - 0.5, color="#0b0b0b", lw=1.0)
ax.set_yticks(range(len(rows_order)))
ax.set_yticklabels([r.replace("Energy/lactate metabolism", "Energy/lactate metab.")
                    .replace("Glutamate-glutamine cycle", "Glu-Gln cycle")
                    for r in rows_order], fontsize=7.5)
ax.set_xticks(range(3))
ax.set_xticklabels(["MERFISH", "Stereo-seq", "Xenium"], fontsize=8)
ax.text(-0.58, (n_func - 1) / 2, "functional", transform=ax.get_yaxis_transform(),
        rotation=90, va="center", fontsize=8, fontweight="bold", color="#52514e")
ax.text(-0.58, n_func + (len(rows_order) - n_func - 1) / 2,
        "metabolic", transform=ax.get_yaxis_transform(),
        rotation=90, va="center", fontsize=8, fontweight="bold", color="#52514e")
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("module Δ (FC − BG)", fontsize=8)
ax.set_title("Functional & metabolic modules — cross-platform summary\n"
             "(* padj < 0.05 within dataset; n.m. = not measurable; full table released)",
             loc="left", fontsize=9, fontweight="bold")
save(fig, "FigureS6_metabolic_matched")
plt.close(fig)

# ================= S8 — per-section lean =================
ps = pd.read_csv(TAB + "/harmonized_lean_persection.csv")
fig, ax = plt.subplots(figsize=(W2, 3.4))
XPOS = {"zhuang": 0, "hao": 1, "suyama": 2}
AN_MARK = {}
markers = ["o", "s", "^", "D", "v"]
rng = np.random.default_rng(1)
for ds in DS_ORDER:
    sub = ps[ps.dataset == ds]
    animals = sorted(sub.animal.unique())
    for ai, an in enumerate(animals):
        aa = sub[sub.animal == an]
        jit = rng.uniform(-0.16, 0.16, len(aa))
        ax.scatter(XPOS[ds] + jit, aa.delta, s=8 + aa.n_pairs * 1.1,
                   marker=markers[ai % len(markers)], facecolors="none",
                   edgecolors=GROUP["FC"] if ds != "hao" else "#52514e", linewidths=0.9,
                   label=None)
    med = np.average(sub.delta, weights=sub.n_pairs)
    ax.hlines(med, XPOS[ds] - 0.28, XPOS[ds] + 0.28, color="k", lw=1.8, zorder=5)

ax.axhline(0, color="#898781", lw=0.8, ls="--")
ax.set_xticks(range(3)); ax.set_xticklabels([DSL[k] for k in DS_ORDER], fontsize=8)
ax.set_ylabel("per-section paired Δ astro module\n(FC − BG)")
ax.set_ylim(-0.78, 0.85)

save(fig, "FigureS8_persection_lean")
plt.close(fig)

# S10 (selection-margin / manual-review sensitivity) is built by FigureS10_threshold.py from
# tables/selection_sensitivity.csv (build_selection_sensitivity.py) since the November 2026 review.
