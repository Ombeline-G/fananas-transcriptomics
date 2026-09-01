#!/usr/bin/env python
"""Figure 5 — 'A modest astrocytic expression shift accompanies the ectopic soma position'
(retitled per review point 13; 'replication' is reserved for Aqp4/Gfap on the platforms
measuring them comparably — decision A).
a: the shift at its three levels of structure (review point 16) — individual matched
   pairs in the background, section means, animal means as the principal level;
   no section-level p in this figure (sections are nested within animals).
b: per-gene matched delta (forest): Aqp4/S100b/Gfap carry the lean; Agt discordant; Hao gate genes flagged.
c: parenchymal-astrocyte PROXIMITY control, per section (review points 20, 25 —
   excludes one specific source only): parenchymal astrocytes sit closer to BG
   than to FC in nearly every section (open markers) and in the pooled medians
   (filled).
All values are read from the canonical tables (harmonized_lean_*.csv, fig4S3_spillover.csv).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, DSL, DS_ORDER, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"

C_NS = "#b8b6b0"

fig, axes = plt.subplots(1, 3, figsize=(W2, 3.35), gridspec_kw=dict(wspace=0.44))
fig.subplots_adjust(top=0.84, bottom=0.25)
axA, axB, axC = axes

# ---------- a: module deltas — cells, sections, animals (review point 16) ----------
mod = pd.read_csv(TAB + "/harmonized_lean_module.csv").set_index("dataset").loc[DS_ORDER]
persec = pd.read_csv(TAB + "/harmonized_lean_persection.csv")
pairs = pd.read_csv(TAB + "/harmonized_lean_pairs.csv")
assert mod.n_pairs.astype(int).tolist() == [335, 46, 42]   # 4 animals (ABCA-1 integrated 2026-08-18)
assert np.allclose(mod.delta_module, [0.17024300, 0.01162437, 0.24434682], atol=1e-4)
assert pairs.groupby("dataset").size().to_dict() == {"zhuang": 335, "hao": 46, "suyama": 42}
x = np.arange(3)
rng = np.random.default_rng(0)
for xi, ds in enumerate(DS_ORDER):
    pp = pairs[pairs.dataset == ds]
    axA.scatter(xi + rng.uniform(-0.24, 0.16, len(pp)), pp.pair_diff, s=4, c="#dddbd4",
                linewidths=0, zorder=2, label="matched pair" if xi == 0 else None)
    ss = persec[persec.dataset == ds]
    axA.scatter(xi + rng.uniform(-0.20, 0.12, len(ss)), ss.delta, s=15, facecolors="none",
                edgecolors="#52514e", linewidths=0.8, zorder=4,
                label="section mean" if xi == 0 else None)
    aa = pp.groupby("animal").pair_diff.mean().sort_values()
    axA.scatter(np.linspace(0.26, 0.36, len(aa)) + xi, aa.values, marker="D", s=34,
                c=GROUP["FC"], edgecolors="#333333", linewidths=0.6, zorder=6,
                label="animal mean" if xi == 0 else None)
za = pairs[pairs.dataset == "zhuang"].groupby("animal").pair_diff.mean()
assert np.allclose(sorted(za.round(2)), [0.12, 0.12, 0.18, 0.22], atol=0.005)   # 4/4 animals positive
axA.axhline(0, color="k", lw=0.8)
axA.set_xticks(x); axA.set_xticklabels(["MERFISH", "Stereo-seq", "Xenium"], fontsize=7, rotation=90, va="top")
axA.set_ylabel("astrocyte-module Δ (FC − BG)\ncapture-matched pairs")
axA.set_ylim(-1.55, 2.3)
axA.legend(frameon=False, fontsize=7, loc="upper right", labelspacing=0.3,
           handletextpad=0.25, borderaxespad=0.02, markerscale=1.2)
axA.set_title("Astrocyte-module lean", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axA, "a", dx=-0.22, dy=1.12)

# ---------- b: per-gene forest ----------
genes = ["Aqp4", "S100b", "Gfap", "Ntsr2", "Aldh1l1", "Gja1", "Agt"]
gtab = pd.read_csv(TAB + "/harmonized_lean_genes.csv")
DCOL = {"zhuang": "#5e3c99", "hao": "#898781", "suyama": "#17becf"}
DMRK = {"zhuang": "o", "hao": "^", "suyama": "s"}
OFF = {"zhuang": -0.22, "hao": 0.0, "suyama": 0.22}
for gi, g in enumerate(genes):
    for ds in DS_ORDER:
        r = gtab[(gtab.dataset == ds) & (gtab.gene == g)]
        if r.empty:
            continue
        r = r.iloc[0]
        gate = r.interpretable.startswith("NO")
        sig = (r.p_wilcoxon < 0.05) and not gate
        style = dict(marker="x", color=DCOL[ds], s=30) if gate else dict(
            marker=DMRK[ds], facecolors=DCOL[ds] if sig else "none",
            edgecolors=DCOL[ds], s=26)
        axB.scatter(r.delta_z, gi + OFF[ds], linewidths=1.0, zorder=4, **style)
axB.axvline(0, color="k", lw=0.8)
axB.set_yticks(range(len(genes)))
axB.set_yticklabels([f"$\\it{{{g}}}$" for g in genes], fontsize=8)
axB.set_ylim(len(genes) - 0.5, -0.9)
axB.set_xlabel("per-gene Δz (FC − BG), matched pairs")
axB.set_xlim(-0.45, 0.95)
h = [plt.Line2D([0], [0], marker=DMRK[d], color="w", markerfacecolor=DCOL[d], markersize=5.5,
                label={"zhuang": "MERFISH", "hao": "Stereo-seq", "suyama": "Xenium"}[d]) for d in DS_ORDER]
h += [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#52514e",
                 markersize=5.5, label="open = n.s."),
      plt.Line2D([0], [0], marker="x", color="#898781", lw=0, markersize=5.5, label="× = gate gene (Stereo-seq)")]
fig.legend(handles=h, frameon=False, fontsize=7, loc="lower center", ncol=5,
           handletextpad=0.3, columnspacing=0.9, bbox_to_anchor=(0.5, 0.015))
axB.set_title("Carrier genes", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axB, "b", dx=-0.16, dy=1.12)

# ---------- c: parenchymal proximity, per section (review points 20, 25) ----------
sp = pd.read_csv(TAB + "/fig4S3_spillover.csv")
sp["key"] = ["zhuang", "hao", "suyama"]
sp = sp.set_index("key").loc[DS_ORDER]
assert sp.n_FC.astype(int).tolist() == [358, 89, 44]
assert np.allclose(sp.dpar_FC_um, [302.8284, 276.5863, 71.6495], atol=1e-3)
assert np.allclose(sp.dpar_BG_um, [243.9175, 210.2380, 38.7451], atol=1e-3)
# 4 animals: within-FC rho(score, d_parench) is now weakly significant in MERFISH
# (rho -0.111, p 0.035, n=358, cell-level descriptive) — declared in text; the first-order
# geometric control is unchanged (parenchymal astrocytes lie FARTHER from FC than from BG).
assert (sp.p_prox.iloc[1:] > 0.05).all() and abs(float(sp.spearman_astro_vs_dist.iloc[0]) + 0.111) < 0.01
persec = pd.read_csv(TAB + "/spillover_persection.csv")
MRK = {"zhuang": "o", "hao": "^", "suyama": "s"}
LBL = {"zhuang": "MERFISH section", "hao": "Stereo-seq section", "suyama": "Xenium (1 section)"}
lims = [8, 2000]
axC.plot(lims, lims, ls="--", c="#898781", lw=0.9)
for ds in DS_ORDER:
    pp = persec[persec.dataset == ds]
    axC.scatter(pp.median_BG_um, pp.median_FC_um, marker=MRK[ds], s=22,
                facecolors="none", edgecolors="#52514e", linewidths=0.9, zorder=4,
                label=LBL[ds])
    axC.scatter([sp.loc[ds, "dpar_BG_um"]], [sp.loc[ds, "dpar_FC_um"]], marker=MRK[ds],
                s=42, c=GROUP["FC"], edgecolors="#333333", linewidths=0.6, zorder=6)
axC.set_xscale("log"); axC.set_yscale("log")
axC.set_xlim(lims); axC.set_ylim(lims)
above = int((persec.median_FC_um > persec.median_BG_um).sum())
axC.text(0.03, 0.97, f"above line = parenchymal\nastro farther from FC\n({above}/{len(persec)} sections;\nfilled = pooled median)",
         zorder=9, path_effects=[pe.withStroke(linewidth=2.4, foreground="white")],
         transform=axC.transAxes, va="top", fontsize=7, color="#52514e")
axC.set_xlabel("BG: median distance to nearest\nparenchymal astrocyte (µm)")
axC.set_ylabel("FC: median distance (µm)")
axC.legend(frameon=False, fontsize=7, loc="lower right", labelspacing=0.3,
           handletextpad=0.2, borderaxespad=0.05)
axC.set_title("Parenchymal proximity control", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axC, "c", dx=-0.22, dy=1.12)

save(fig, "Figure5_lean")
