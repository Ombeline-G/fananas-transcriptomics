#!/usr/bin/env python
"""Figure S12 — dissociated-atlas consistency check (demoted from Figure 4,
review point 12: supplementary support, not a pillar).

The unsupervised part is score-free (Leiden reclustering of Bergmann nuclei);
the FC-shift score is used only to interrogate the resulting subclusters, and
can by construction detect only a state carrying the astrocytic shift observed
in situ — the one measurable FC-distinguishing profile. A hypothetical
dissociated FC state distinguished by other genes would escape this test.

a: per-subcluster FC-shift scores (every Bergmann subcluster of each atlas)
   against the astrocyte reference mean.
b: robustness to batch structure: max subcluster lean recomputed within each
   human donor and each mouse sex, versus the astrocyte mean of the same
   stratum.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT
TAB = str(PROJECT_ROOT) + "/comparison/tables"

fig, (axA, axB) = plt.subplots(1, 2, figsize=(W2 * 0.92, 3.0),
                               gridspec_kw=dict(wspace=0.45, top=0.86, bottom=0.20,
                                                width_ratios=[0.82, 1.18]))

sub_h = pd.read_csv(TAB + "/atlas_siletti_bergmann_subclusters.csv")
sub_m = pd.read_csv(TAB + "/atlas_kozareva_bergmann_subclusters.csv")
assert len(sub_h) == 6 and len(sub_m) == 3
assert np.allclose([sub_h.S_FC_astro_lean.max(), sub_m.S_FC_astro_lean.max()],
                   [-0.29490045, -0.36326173], atol=1e-7)
ASTRO_REF = {"siletti": 0.66, "kozareva": 0.32}  # CANONICAL_NUMBERS section 4

# ---- a: every Bergmann subcluster vs the astrocyte reference ----
rng = np.random.default_rng(0)
for xpos, (key, sub, lab) in enumerate([("siletti", sub_h, "Siletti\n(human)"),
                                        ("kozareva", sub_m, "Kozareva\n(mouse)")]):
    jit = rng.uniform(-0.10, 0.10, len(sub))
    axA.scatter(xpos + jit, sub.S_FC_astro_lean, s=26, facecolors="none",
                edgecolors="#52514e", linewidths=1.0, zorder=4,
                label="Bergmann subclusters" if xpos == 0 else None)
    axA.hlines(ASTRO_REF[key], xpos - 0.22, xpos + 0.22, color=GROUP["OtherCBastro"],
               lw=2.2, zorder=5, label="astrocyte reference" if xpos == 0 else None)
axA.axhline(0, color="k", lw=0.8)
axA.set_xticks([0, 1]); axA.set_xticklabels(["Siletti\n(human)", "Kozareva\n(mouse)"], fontsize=8)
axA.set_xlim(-0.55, 1.55)
axA.set_ylabel("FC-shift score (z)")
axA.set_ylim(-0.75, 0.95)
axA.legend(frameon=False, fontsize=7, loc="upper left", labelspacing=0.3)
axA.set_title("All Bergmann subclusters", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axA, "a", dx=-0.24, dy=1.10)

# ---- b: robustness within donors (human) and sexes (mouse) ----
don = pd.read_csv(TAB + "/atlas_siletti_donor_verdict.csv")
sex = pd.read_csv(TAB + "/atlas_kozareva_sex_verdict.csv")
assert len(don) == 3 and len(sex) == 2
strata = ([(f"donor {i + 1}", r.max_subcluster_lean, r.astro_mean_lean)
           for i, r in don.iterrows()]
          + [(f"{r.sex}", r.max_subcluster_lean, r.astro_mean_lean) for _, r in sex.iterrows()])
xs = np.arange(len(strata)); wb = 0.34
axB.bar(xs - wb / 2, [s[1] for s in strata], wb, color="#52514e",
        label="max subcluster lean")
axB.bar(xs + wb / 2, [s[2] for s in strata], wb, color=GROUP["OtherCBastro"],
        label="astrocyte mean")
axB.axhline(0, color="k", lw=0.8)
axB.axvline(2.5, color="#c3c2b7", lw=0.8)
axB.text(1.0, 0.92, "human donors", ha="center", fontsize=7.5, color="#52514e")
axB.text(3.5, 0.92, "mouse sexes", ha="center", fontsize=7.5, color="#52514e")
axB.set_xticks(xs)
axB.set_xticklabels([s[0] for s in strata], fontsize=7.5, rotation=35, ha="right")
axB.set_ylabel("FC-shift score (z)")
axB.set_ylim(-0.75, 1.05)
axB.legend(frameon=False, fontsize=7, loc="lower right", labelspacing=0.3)
axB.set_title("Within batch strata", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axB, "b", dx=-0.24, dy=1.10)

save(fig, "FigureS12_atlases")
