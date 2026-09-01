#!/usr/bin/env python
"""Figure S15 — vascular proximity does not explain the astrocytic accent
(review point 21; negative for the perivascular hypothesis, reported with its tests).

Soma-to-vessel distance DOES differ significantly from BG on both single-soma
platforms — but in OPPOSITE directions (MERFISH: FC closer, 36 vs 43 um,
section-level paired p = 0.019; Xenium: FC farther, 32 vs 25 um, p = 2.8e-4).
The platforms therefore do not support a consistent perivascular FC position,
and within FC Aqp4 does not track vessel distance on either platform.

a: MERFISH — per-section median distance to the nearest endothelial soma,
   FC versus BG (identity line; FC closer in 24/34 sections, paired p = 0.019).
b: Xenium — per-cell distances to the Vascular/stroma proxy, FC versus BG
   (single section; FC farther, cell-level p = 2.8e-4, descriptive).
c: MERFISH, within FC — Aqp4 (z) against vessel distance (descriptive rho).
Distances recomputed deterministically here and asserted against
vascular_proximity.csv.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT, ZHUANG, SUYAMA

TAB = str(PROJECT_ROOT) + "/comparison/tables"
canon = pd.read_csv(TAB + "/vascular_proximity.csv")

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(W2, 3.5),
                                    gridspec_kw=dict(wspace=0.52, top=0.88, bottom=0.28))

# ---- recompute MERFISH distances ----
m = pd.read_parquet(TAB + "/zhuang_astro_cells.parquet")
E = pd.read_parquet(TAB + "/zhuang_astro_expr.parquet").reindex(m.cell)
z = (E - E.mean()) / E.std(ddof=0)
keep = m.group.isin(["FC", "BG"]).values
b = m[keep].copy()
b["aqp4"] = z["Aqp4"].values[keep]
cb = pd.read_parquet(str(ZHUANG / "cb_cells_layers_all4.parquet"))
endo = cb[cb.subclass == "333 Endo NN"]
dists = np.full(len(b), np.nan)
for sec, sub in b.groupby("section"):
    ev = endo[endo.brain_section_label == sec]
    tree = cKDTree(ev[["x", "y"]].to_numpy())
    d, _ = tree.query(sub[["x", "y"]].to_numpy())
    dists[b.index.get_indexer(sub.index)] = d * 1000.0
b["dvasc"] = dists
summ = canon[(canon.dataset == "zhuang") & (canon.level == "summary")].iloc[0]
assert abs(b[b.group == "FC"].dvasc.median() - summ.median_FC_um) < 0.01
assert abs(b[b.group == "BG"].dvasc.median() - summ.median_BG_um) < 0.01

# ---- a: per-section medians ----
sec = canon[(canon.dataset == "zhuang") & (canon.level == "section")]
assert len(sec) == 34 and int((sec.median_FC_um < sec.median_BG_um).sum()) == 24   # 4 animals
lim = [0, 105]
axA.plot(lim, lim, ls="--", c="#898781", lw=0.9)
axA.scatter(sec.median_BG_um, sec.median_FC_um, s=22, facecolors="none",
            edgecolors="#52514e", linewidths=1.0, zorder=4)
axA.set_xlim(lim); axA.set_ylim(lim)
axA.set_xlabel("BG: median distance to nearest\nendothelial soma (µm)")
axA.set_ylabel("FC: median distance (µm)")
p_sec = float(summ.p_section_paired)
axA.text(0.04, 0.96, f"below line = FC closer\n(24/34 sections;\npaired p = {p_sec:.3f})",
         transform=axA.transAxes, va="top", fontsize=7.2, color="#52514e")
axA.set_title("MERFISH — per section", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axA, "a", dx=-0.24, dy=1.12)

# ---- b: Xenium per-cell strips ----
sm = pd.read_parquet(TAB + "/suyama_astro_cells.parquet")
if "cell" not in sm.columns:
    sm = sm.rename(columns={"index": "cell"})
sm["cell"] = sm.cell.astype(str)
sb = sm[sm.group.isin(["FC", "BG"])].copy()
sc = pd.read_csv(str(SUYAMA / "results" / "tables" / "adult_cells_scores.csv"), dtype={"cell": str})
vasc = sc[sc.cl.isin(["c8", "c20", "c21", "c23"])]
tree = cKDTree(vasc[["x", "y"]].to_numpy())
sb["dvasc"], _ = tree.query(sb[["x", "y"]].to_numpy())
ssum = canon[(canon.dataset == "suyama") & (canon.level == "summary")].iloc[0]
assert abs(sb[sb.group == "FC"].dvasc.median() - ssum.median_FC_um) < 0.01
rng = np.random.default_rng(0)
for xi, grp in enumerate(["BG", "FC"]):
    v = sb[sb.group == grp].dvasc
    axB.scatter(xi + rng.uniform(-0.18, 0.18, len(v)), v, s=4 if grp == "BG" else 10,
                c=GROUP[grp], alpha=0.4, linewidths=0)
    axB.hlines(v.median(), xi - 0.26, xi + 0.26, color="#0b0b0b", lw=1.6, zorder=5)
axB.set_xticks([0, 1]); axB.set_xticklabels(["BG", "FC"], fontsize=8)
axB.set_xlim(-0.6, 1.6)
axB.set_ylabel("distance to Vascular/stroma\nproxy (µm)")
axB.text(0.97, 0.97, f"FC farther\np = {ssum.p_MW_cell:.0e} (cell-level,\ndescriptive)",
         transform=axB.transAxes, ha="right", va="top", fontsize=7.2, color="#52514e")
axB.set_title("Xenium — single section", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axB, "b", dx=-0.26, dy=1.12)

# ---- c: within-FC Aqp4 vs vessel distance (MERFISH) ----
fc = b[b.group == "FC"]
rho = spearmanr(fc.aqp4, fc.dvasc)[0]
assert abs(rho - summ.rho_aqp4_dvasc_FC) < 1e-6
axC.scatter(fc.dvasc, fc.aqp4, s=10, c=GROUP["FC"], alpha=0.55, linewidths=0)
axC.set_xlabel("FC distance to nearest\nendothelial soma (µm)")
axC.set_ylabel("Aqp4 (z)")
axC.text(0.96, 0.96, f"ρ = {rho:+.2f} (n = {len(fc)};\ndescriptive)",
         transform=axC.transAxes, ha="right", va="top", fontsize=7.2, color="#52514e")
axC.set_title("MERFISH — within FC", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axC, "c", dx=-0.24, dy=1.12)

fig.text(0.5, 0.035, "Significant on both platforms but in opposite directions — no consistent perivascular position;\n"
         "within FC, Aqp4 does not track vessel distance (c).",
         ha="center", va="top", fontsize=7.2, color="#52514e")
save(fig, "FigureS15_vessels")
