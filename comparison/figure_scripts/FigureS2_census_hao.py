#!/usr/bin/env python
"""Figure S2 — Mouse1 census robustness on the whole-transcriptome Hao platform + historical markers.
Demoted from main Figure 2 (editorial pass 2026-08-10): panel a is a caveat (ML non-BG somata are
low-capture, no coherent identity), panel b a null on the contested historical markers — supporting
material, not main-figure results.
a: identity scores as point-range (median ± bootstrap 95% CI) for the four soma groups;
b: % of cells expressing Kcnip3/Kcnb2 with Wilson 95% CIs.
G1 audit asserts all values against CANONICAL_NUMBERS; aborts on mismatch.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, panel_letter, save, W2
from panels import module_scores

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"

SCORE_PAL = {"berg": "#5e3c99", "astro": "#e08214"}
GENE_COL = {"Kcnip3": "#2166ac", "Kcnb2": "#999999"}
GORD = ["BG", "FC", "Astro_PCL_nonBG", "Astro_ML_nonBG"]
GLAB = ["PCL", "ML", "PCL", "ML"]


def wilson(k, n, z=1.959964):
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    hw = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return 100 * (ctr - hw), 100 * (ctr + hw)


# ================= G1 — DATA AUDIT =================
print("== G1 DATA AUDIT (Figure S2) ==")
hs = pd.read_csv(TAB + "/hao_ML_census_summary.csv").set_index("group")
assert [int(hs.loc[g, "n"]) for g in GORD] == [3723, 450, 414, 114]
assert [round(hs.loc[g, "nCount_med"]) for g in GORD] == [2083, 1562, 1654, 883]
assert [round(hs.loc[g, "Kcnb2_pct"], 1) for g in GORD] == [0.8, 1.3, 0.5, 0.0]
assert round(hs.loc["Astro_ML_nonBG", "Kcnip3_pct"], 1) == 6.1
assert round(hs.loc["FC", "Kcnip3_pct"], 1) == 6.9
he = pd.read_parquet(TAB + "/hao_ML_census_expr.parquet")
assert set(he.mouse.astype(str)) == {"Mouse1"}
sc, _ = module_scores(he.select_dtypes(include="number").drop(
    columns=[c for c in ("nCount",) if c in he.columns]))
sc["group"] = he["group"].values
for g in GORD:
    assert abs(float(sc.loc[sc.group == g, "bergmann"].median()) - float(hs.loc[g, "berg_score"])) < 0.01
    assert abs(float(sc.loc[sc.group == g, "astrocyte"].median()) - float(hs.loc[g, "astro_score"])) < 0.01
print("  Mouse1 only ✓  n 3723/450/414/114 ✓  tx 2083/1562/1654/883 ✓  Kcnip3 6.1/6.9 ✓  scores ✓")

# ================= figure =================
fig = plt.figure(figsize=(W2, 3.3))
gs = GridSpec(1, 2, figure=fig, wspace=0.42, top=0.86, bottom=0.30, left=0.10, right=0.975)
axC = fig.add_subplot(gs[0, 0]); axD = fig.add_subplot(gs[0, 1])
fig.text(0.5, 0.955, "Stereo-seq Mouse1, position-first census", ha="center",
         fontsize=9.5, fontweight="bold")


def brackets(ax, y_rule=-0.175, y_lab=-0.21):
    for x0, x1, lab in [(0, 1, "Bergmann somata"), (2, 3, "non-BG astrocytes")]:
        tr = ax.get_xaxis_transform()
        ax.plot([x0 - 0.30, x1 + 0.30], [y_rule, y_rule], transform=tr, color="#52514e",
                lw=0.8, clip_on=False)
        for xe in (x0 - 0.30, x1 + 0.30):
            ax.plot([xe, xe], [y_rule, y_rule + 0.02], transform=tr, color="#52514e",
                    lw=0.8, clip_on=False)
        ax.text((x0 + x1) / 2, y_lab, lab, transform=tr, ha="center", va="top",
                fontsize=8, fontweight="bold")


# ---------- a: identity scores, point-range ----------
xg = np.arange(4)
OFF = {"berg": -0.13, "astro": 0.13}
KEY = {"berg": "bergmann", "astro": "astrocyte"}
rngb = np.random.default_rng(0)
for score, off in OFF.items():
    med, lo, hi = [], [], []
    for g in GORD:
        v = sc.loc[sc.group == g, KEY[score]].values
        boots = np.median(rngb.choice(v, size=(500, len(v)), replace=True), axis=1)
        med.append(np.median(v)); lo.append(np.percentile(boots, 2.5)); hi.append(np.percentile(boots, 97.5))
    med, lo, hi = map(np.array, (med, lo, hi))
    axC.errorbar(xg + off, med, yerr=[med - lo, hi - med], fmt="o", ms=4.5,
                 color=SCORE_PAL[score], ecolor=SCORE_PAL[score], elinewidth=1.0, capsize=2)
axC.axhline(0, color="#4d4d4d", lw=0.8)
axC.set_xticks(xg); axC.set_xticklabels(GLAB, fontsize=8)
axC.tick_params(axis="x", length=0)
axC.set_xlim(-0.6, 3.6)
axC.set_ylabel("identity score\n(whole-transcriptome z, median)")
axC.set_ylim(-0.68, 0.12); axC.set_yticks([0.0, -0.2, -0.4, -0.6])
brackets(axC)
hC = [Line2D([0], [0], marker="o", color=SCORE_PAL["berg"], lw=0, markersize=5,
             label="Bergmann identity"),
      Line2D([0], [0], marker="o", color=SCORE_PAL["astro"], lw=0, markersize=5,
             label="astrocyte identity")]
axC.legend(handles=hC, frameon=False, fontsize=7, loc="upper center",
           bbox_to_anchor=(0.5, -0.34), ncol=2, columnspacing=1.2, handletextpad=0.35)
panel_letter(axC, "a", dx=-0.20, dy=1.10)

# ---------- b: % expressing with Wilson CIs ----------
w = 0.30
for j, gene in enumerate(["Kcnip3", "Kcnb2"]):
    vals = [float(hs.loc[g, f"{gene}_pct"]) for g in GORD]
    ns = [int(hs.loc[g, "n"]) for g in GORD]
    ks = [round(v * n / 100) for v, n in zip(vals, ns)]
    cis = [wilson(k, n) for k, n in zip(ks, ns)]
    errd = np.array([[v - lo, hi - v] for v, (lo, hi) in zip(vals, cis)]).T
    xb = xg + (j - 0.5) * w
    axD.bar(xb, vals, w, color=GENE_COL[gene], alpha=0.9, edgecolor="#333333", linewidth=0.6)
    axD.errorbar(xb, vals, yerr=errd, fmt="none", ecolor="#0b0b0b", elinewidth=0.8, capsize=2)
axD.set_xticks(xg); axD.set_xticklabels(GLAB, fontsize=8)
axD.tick_params(axis="x", length=0)
axD.set_xlim(-0.6, 3.6)
axD.set_ylabel("cells expressing (%)")
axD.set_ylim(0, 15); axD.set_yticks([0, 5, 10, 15])
brackets(axD)
hD = [plt.Rectangle((0, 0), 1, 1, facecolor=GENE_COL["Kcnip3"], label="Kcnip3 (KChIP3)"),
      plt.Rectangle((0, 0), 1, 1, facecolor=GENE_COL["Kcnb2"], label="Kcnb2 (Kv2.2)")]
axD.legend(handles=hD, frameon=False, fontsize=7, loc="upper center",
           bbox_to_anchor=(0.5, -0.34), ncol=2, columnspacing=1.2, handletextpad=0.35)
panel_letter(axD, "b", dx=-0.16, dy=1.10)

save(fig, "FigureS2_census_hao")
