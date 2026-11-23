#!/usr/bin/env python
"""Figure S14 — astrocytic score versus soma position across all Bergmann-class
cells (review point 19).

a: MERFISH — per-cell score against signed perpendicular distance to the local
   Purkinje line (NN units; negative = GCL side), FC red / BG blue, binned
   means overlaid; the GCL side is flat (negative control), the ML side rises;
   the 1-6-unit zone is sparse by census construction (threshold ~45 um).
b: Xenium — same display on the unsigned soma distance of the single section;
   no depth dependence is detectable there.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT, ZHUANG

TAB = str(PROJECT_ROOT) + "/comparison/tables"
MODULES = {"zhuang": ["Aqp4", "Agt", "Gja1", "Aldh1l1"],
           "suyama": ["Aqp4", "Agt", "Gja1", "Slc1a3", "S100b", "Aldh1l1", "Gjb6"]}


def load(key):
    m = pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    E = pd.read_parquet(f"{TAB}/{key}_astro_expr.parquet").reindex(m.cell)
    z = (E - E.mean()) / E.std(ddof=0)
    keep = m.group.isin(["FC", "BG"]).values
    b = m[keep].copy()
    b["score"] = z.loc[:, MODULES[key]].mean(axis=1).values[keep]
    return b


prof = pd.read_csv(TAB + "/lean_depth_continuum_profile.csv")
zp = prof[(prof.dataset == "zhuang") & (prof.animal == "all")]
assert int(zp.n.sum()) == 6446   # 4 animals
assert abs(float(zp[zp.bin == "PCL-1..1"].mean_score.iloc[0]) - (-0.250)) < 0.005
assert abs(float(zp[zp.bin == ">15"].mean_score.iloc[0]) - (-0.138)) < 0.005

fig, (axA, axB) = plt.subplots(1, 2, figsize=(W2, 3.4),
                               gridspec_kw=dict(width_ratios=[1.5, 1], wspace=0.42,
                                                top=0.86, bottom=0.16))

# ---- a: MERFISH ----
b = load("zhuang").set_index("cell")
cb = pd.read_parquet(str(ZHUANG / "cb_cells_layers_all4.parquet"))
b = b.join(cb[["mli_frac", "gran_frac"]])
b["sdepth"] = np.where(b.mli_frac > b.gran_frac, 1.0, -1.0) * b.perp_norm
x = b.sdepth.clip(-20, 40)
for grp, size, zo in [("BG", 3, 2), ("FC", 10, 3)]:
    mm = b.group == grp
    axA.scatter(x[mm], b.score[mm], s=size, c=GROUP[grp], alpha=0.35, linewidths=0, zorder=zo)
BINS = [-np.inf, -1, 1, 3, 6, 10, 15, np.inf]
CENTERS = [-6, 0, 2.4, 4.5, 8, 12.5, 20]
bin_idx = pd.cut(b.sdepth, BINS, labels=False)
means = b.groupby(bin_idx).score.mean()
ns = b.groupby(bin_idx).score.size()
axA.plot([CENTERS[int(i)] for i in means.index], means.values, "-", color="#0b0b0b",
         lw=1.4, zorder=5)
axA.scatter([CENTERS[int(i)] for i in means.index], means.values, s=30, c="#0b0b0b",
            zorder=6, label="bin mean")
for i in means.index:
    yrow = 1.66 if int(i) % 2 == 0 else 1.42
    axA.text(CENTERS[int(i)], yrow, str(int(ns[i])), ha="center",
             fontsize=7, color="#52514e")
axA.text(-19.5, 1.66, "n per bin:", ha="left", fontsize=7, color="#52514e")
axA.axvline(0, color="#b8560f", lw=1.0)
axA.axvspan(1, 6, color="#e1e0d9", alpha=0.45, lw=0)
axA.text(3.5, -2.15, "sparse by census\nconstruction", ha="center", fontsize=7, color="#52514e")
axA.text(-10, -2.15, "GCL side", ha="center", fontsize=7.5, color="#52514e")
axA.text(25, -2.15, "ML side", ha="center", fontsize=7.5, color="#52514e")
axA.set_xlabel("signed distance to local Purkinje line\n(local nearest-neighbor units; clipped at −20/+40)")
axA.set_ylabel("astro-module score (z)")
axA.set_xlim(-21, 41); axA.set_ylim(-2.45, 1.85)
axA.set_title("MERFISH — all Bergmann-class cells", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axA, "a", dx=-0.14, dy=1.12)

# ---- b: Xenium ----
s = load("suyama")
for grp, size, zo in [("BG", 3, 2), ("FC", 10, 3)]:
    mm = s.group == grp
    axB.scatter(s.depth_um[mm], s.score[mm], s=size, c=GROUP[grp], alpha=0.35,
                linewidths=0, zorder=zo, label=grp)
XB = [0, 15, 30, 45, 70, 140]
XC = [7.5, 22.5, 37.5, 57.5, 100]
xb = pd.cut(s.depth_um, XB, labels=False)
xm = s.groupby(xb).score.mean()
axB.plot([XC[int(i)] for i in xm.index], xm.values, "-", color="#0b0b0b", lw=1.4, zorder=5)
axB.scatter([XC[int(i)] for i in xm.index], xm.values, s=30, c="#0b0b0b", zorder=6)
axB.set_xlabel("soma distance to PCL (µm)")
axB.set_ylabel("astro-module score (z)")
axB.set_xlim(0, 140); axB.set_ylim(-2.45, 1.85)
axB.legend(frameon=False, fontsize=7, loc="upper right", labelspacing=0.3, markerscale=1.6)
axB.set_title("Xenium — single section", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axB, "b", dx=-0.22, dy=1.12)

save(fig, "FigureS14_continuum")
