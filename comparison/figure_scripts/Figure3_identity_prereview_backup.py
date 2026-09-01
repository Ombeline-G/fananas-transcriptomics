#!/usr/bin/env python
"""Figure 3 — FC retain Bergmann identity.
a: deterministic kNN — % of each FC's 15 nearest astroglial neighbours that are canonical BG vs the BG
   base rate (fig2_knn_fc_in_bergmann.csv; exactly reproducible, no stochastic step).
b: bootstrap centroid-distance ratios d(FC,other)/d(FC,BG) with 95% CI (fig_bootstrap_ratios.csv; >1 = FC closer to BG).
c: astroglia UMAP per dataset (display only; fixed seed; the claim is panel a, stated in the caption).
d: shared Bergmann signature dotplot (computed from the canonical parquets; grey 'n.m.' = not measured).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, GROUP_LABEL, DSL, DS_ORDER, W2, panel_letter, save, seq_blue

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"

fig = plt.figure(figsize=(W2, 8.6))
gs = GridSpec(3, 6, figure=fig, height_ratios=[0.85, 1.0, 1.05], hspace=0.52, wspace=1.1,
              top=0.95, bottom=0.06)
axA = fig.add_subplot(gs[0, 0:3])
axB = fig.add_subplot(gs[0, 3:6])
axU = [fig.add_subplot(gs[1, 2 * i:2 * i + 2]) for i in range(3)]
axD = fig.add_subplot(gs[2, 0:6])

# ---------- a: kNN ----------
knn = pd.read_csv(TAB + "/fig2_knn_fc_in_bergmann.csv").set_index("key").loc[DS_ORDER]
assert knn.n_FC.astype(int).tolist() == [190, 89, 44]
assert np.allclose(knn.frac_BG, [98.9122807, 68.0898876, 88.4848485], atol=1e-6)
assert np.allclose(knn.base_BG, [71.4428657, 61.3526570, 50.2480706], atol=1e-6)
x = np.arange(3); w = 0.36
axA.bar(x - w / 2, knn.frac_BG, w, color=GROUP["BG"], label="observed: neighbours that are BG")
axA.bar(x + w / 2, knn.base_BG, w, color="#b8b6b0", label="chance: BG share of astroglia")

axA.set_xticks(x); axA.set_xticklabels([DSL[k].replace(" (", "\n(") for k in DS_ORDER], fontsize=8)
axA.set_ylabel("% of FC's 15 nearest\nastroglial neighbours")
axA.set_ylim(0, 132)
axA.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0.01, 1.00))
panel_letter(axA, "a", dx=-0.20, dy=1.14)

# ---------- b: bootstrap centroid ratios ----------
boot = pd.read_csv(TAB + "/fig_bootstrap_ratios.csv").set_index("dataset").loc[DS_ORDER]
assert np.allclose(boot.med, [2.4627516, 1.1318223, 1.8953412], atol=1e-6)
assert (boot.lo > 1).all() and (boot.hi > boot.med).all()
y = np.arange(3)[::-1]
axB.errorbar(boot.med, y, xerr=[boot.med - boot.lo, boot.hi - boot.med], fmt="o",
             color="#0b0b0b", capsize=3, ms=5, lw=1.2)

axB.axvline(1.0, color="#898781", ls="--", lw=0.9)

axB.set_yticks(y); axB.set_yticklabels(["Zhuang", "Hao", "Suyama"], fontsize=8)
axB.set_xlabel("centroid ratio d(FC, other astro) / d(FC, BG)\n(bootstrap median, 95% CI)")
axB.set_xlim(0.8, 2.85); axB.set_ylim(-0.85, 2.6)
panel_letter(axB, "b", dx=-0.32, dy=1.14)

# ---------- c: UMAPs (display projection, fixed seed) ----------
try:
    from umap import UMAP
    HAVE_UMAP = True
except Exception:
    HAVE_UMAP = False
for i, key in enumerate(DS_ORDER):
    ax = axU[i]
    m = pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    if "cell" not in m.columns:
        m = m.rename(columns={"index": "cell"})
        if "cell" not in m.columns:
            m["cell"] = m.index
    E = pd.read_parquet(f"{TAB}/{key}_astro_expr.parquet")
    E = E.reindex(m["cell"].values) if m["cell"].isin(E.index).all() else E
    expected_groups = {
        "zhuang": {"FC": 190, "BG": 3570, "OtherCBastro": 1427},
        "hao": {"FC": 89, "BG": 1524, "OtherCBastro": 960},
        "suyama": {"FC": 44, "BG": 1823, "OtherCBastro": 1805},
    }
    assert m.group.value_counts().to_dict() == expected_groups[key]
    assert len(E) == len(m)
    Z = ((E - E.mean()) / (E.std() + 1e-9)).fillna(0).values
    pcs = PCA(n_components=min(30, Z.shape[1] - 1), random_state=0).fit_transform(Z)
    emb = UMAP(n_neighbors=60, min_dist=0.7, random_state=0).fit_transform(pcs) if HAVE_UMAP else pcs[:, :2]
    g = m.group.values
    for grp, s, z, mk, ec, lw in [("OtherCBastro", 3, 2, "o", "none", 0), ("BG", 3, 3, "o", "none", 0),
                                  ("FC", 26, 5, "*", "#222222", 0.3)]:
        mm = g == grp
        ax.scatter(emb[mm, 0], emb[mm, 1], s=s, marker=mk, c=GROUP[grp], edgecolors=ec, linewidths=lw)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color("#c3c2b7")
    ax.set_xlabel("UMAP1", fontsize=7); ax.set_ylabel("UMAP2", fontsize=7)
    ax.set_title(DSL[key], fontsize=9, fontweight="bold", pad=3)
    if i == 0:
        panel_letter(ax, "c", dx=-0.12, dy=1.16)

# ---------- d: shared Bergmann signature dotplot (from canonical parquets) ----------
GENES = ["Hopx", "Gpr37l1", "A2m", "Tnc", "Agt", "Gdf10", "Slc1a3", "Fabp7"]
GRPS = ["FC", "BG", "OtherCBastro"]
GLAB = {"FC": "FC", "BG": "BG", "OtherCBastro": "other"}
rows = []
for key in DS_ORDER:
    m = pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    if "cell" not in m.columns:
        m = m.rename(columns={"index": "cell"})
        if "cell" not in m.columns:
            m["cell"] = m.index
    E = pd.read_parquet(f"{TAB}/{key}_astro_expr.parquet")
    E = E.reindex(m["cell"].values) if m["cell"].isin(E.index).all() else E
    for gname in GENES:
        if gname not in E.columns:
            for grp in GRPS:
                rows.append(dict(dataset=key, gene=gname, group=grp, mean=np.nan, frac=np.nan))
            continue
        v = E[gname].values
        for grp in GRPS:
            mm = (m.group == grp).values
            rows.append(dict(dataset=key, gene=gname, group=grp,
                             mean=float(np.nanmean(v[mm])), frac=float((v[mm] > 0).mean())))
dp = pd.DataFrame(rows)
# scale colour 0..1 per gene x dataset (across the 3 groups)
dp["scaled"] = dp.groupby(["dataset", "gene"])["mean"].transform(lambda s: s / (s.max() + 1e-9))
ncol = len(DS_ORDER) * len(GRPS)
for gi, gname in enumerate(GENES):
    for di, key in enumerate(DS_ORDER):
        for pi, grp in enumerate(GRPS):
            r = dp[(dp.dataset == key) & (dp.gene == gname) & (dp.group == grp)].iloc[0]
            xpos = di * 3 + pi
            if np.isnan(r["mean"]):
                if pi == 1:
                    axD.text(xpos, gi, "n.m.", ha="center", va="center", fontsize=7, color="#b8b6b0")
                continue
            axD.scatter(xpos, gi, s=8 + 130 * r.frac, c=[seq_blue(r.scaled)],
                        edgecolors="#52514e", linewidths=0.4)
axD.set_yticks(range(len(GENES)))
axD.set_yticklabels([f"$\\it{{{g}}}$" for g in GENES], fontsize=8)
axD.set_xticks(range(ncol))
axD.set_xticklabels([GLAB[g] for _ in DS_ORDER for g in GRPS], fontsize=7.5)
for di, key in enumerate(DS_ORDER):
    axD.text(di * 3 + 1, -1.15, DSL[key], ha="center", fontsize=8.5, fontweight="bold")
    if di:
        axD.axvline(di * 3 - 0.5, color="#e1e0d9", lw=0.8)
axD.set_ylim(len(GENES) - 0.4, -1.75); axD.set_xlim(-0.6, ncol - 0.4)
# size + colour keys
for j, fr in enumerate([0.25, 0.5, 1.0]):
    axD.scatter(ncol + 0.55, len(GENES) - 3.1 - j * 1.05, s=8 + 130 * fr, c="#b8b6b0",
                edgecolors="#52514e", linewidths=0.4, clip_on=False)
    axD.text(ncol + 0.95, len(GENES) - 3.1 - j * 1.05, f"{int(fr*100)}%", fontsize=7, va="center", clip_on=False)
axD.text(ncol + 0.55, len(GENES) - 2.05, "% expr.", fontsize=7, ha="center", clip_on=False)
for j, sc in enumerate([1.0, 0.5, 0.15]):
    axD.scatter(ncol + 0.55, 1.0 + j * 0.75, s=60, c=[seq_blue(sc)], edgecolors="#52514e",
                linewidths=0.4, clip_on=False)
axD.text(ncol + 0.55, 0.25, "mean expr.\n(scaled)", fontsize=7, ha="center", clip_on=False)
panel_letter(axD, "d", dx=-0.065, dy=1.14)

save(fig, "Figure3_identity")
