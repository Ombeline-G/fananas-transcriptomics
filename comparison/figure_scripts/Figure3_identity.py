#!/usr/bin/env python
"""Figure 3 — FC retain Bergmann identity (reordered 2026-08-17, review point 5).

a: the Bergmann signature is DEFINED WITHOUT FC — dotplot of signature genes in BG vs
   parenchymal astrocytes only (computed from the canonical parquets; 'n.m.' = not measured).
b: FC independently reproduce the signature — FC dotplot on the same genes and colour scale,
   plus Δlog2(BG vs parenchymal) against Δlog2(FC vs parenchymal) per gene x dataset
   (fig6_all_contrasts_DE.csv), with Pearson r.
c: quantitative identity — deterministic kNN (fig2_knn_fc_in_bergmann.csv) and centroid
   ratios d(FC, parenchymal)/d(FC, BG) (>1 = FC closer to BG): per-section replication
   (centroid_ratio_by_section.csv, open circles; conservative — attenuated toward 1 by
   few-FC noise) with the pooled all-cell ratio as a tick; the single-section Xenium
   dataset shows the pooled ratio with its DESCRIPTIVE cell-bootstrap 95% interval
   (fig_bootstrap_ratios.csv — cells are not independent, review point 6).
d: astroglia UMAP per dataset (display only; fixed seed; the claims are panels a-c).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, DSL, DS_ORDER, W2, panel_letter, save, seq_blue

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"

GENES = ["Hopx", "Gpr37l1", "A2m", "Tnc", "Agt", "Gdf10", "Slc1a3", "Fabp7"]
GRPS = ["BG", "OtherCBastro", "FC"]

# ---------- dotplot data (all three groups at once, shared scaling) ----------
rows = []
for key in DS_ORDER:
    m = pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    if "cell" not in m.columns:
        m = m.rename(columns={"index": "cell"})
        if "cell" not in m.columns:
            m["cell"] = m.index
    E = pd.read_parquet(f"{TAB}/{key}_astro_expr.parquet")
    E = E.reindex(m["cell"].values) if m["cell"].isin(E.index).all() else E
    expected = {"zhuang": {"FC": 362, "BG": 6084, "OtherCBastro": 1910},   # 4 animals since ABCA-1 integration (2026-08-18)
                "hao": {"FC": 89, "BG": 1524, "OtherCBastro": 960},
                "suyama": {"FC": 44, "BG": 1823, "OtherCBastro": 1805}}
    assert m.group.value_counts().to_dict() == expected[key]
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
dp["scaled"] = dp.groupby(["dataset", "gene"])["mean"].transform(lambda s: s / (s.max() + 1e-9))

# ---------- contrasts for the reproduction scatter ----------
de = pd.read_csv(TAB + "/fig6_all_contrasts_DE.csv")
de = de[de.gene.isin(GENES)]
bgc = de[de.contrast == "BG_vs_OtherAstro"].set_index(["gene", "dataset"])
fcc = de[de.contrast == "FC_vs_OtherAstro"].set_index(["gene", "dataset"])
common = bgc.index.intersection(fcc.index)
xv = bgc.loc[common, "log2FC"].values; yv = fcc.loc[common, "log2FC"].values
ok = np.isfinite(xv) & np.isfinite(yv)
r = float(np.corrcoef(xv[ok], yv[ok])[0, 1])
print(f"signature reproduction: {ok.sum()} gene x dataset pairs, Pearson r = {r:.3f}")
assert r > 0.8, r

# ================= figure =================
fig = plt.figure(figsize=(W2, 9.2))
gs = GridSpec(3, 6, figure=fig, height_ratios=[1.15, 0.82, 0.95], hspace=0.52, wspace=1.1,
              top=0.955, bottom=0.055)
axA = fig.add_subplot(gs[0, 0:3])
axB1 = fig.add_subplot(gs[0, 3:4])
axB2 = fig.add_subplot(gs[0, 4:6])
axC1 = fig.add_subplot(gs[1, 0:3])
axC2 = fig.add_subplot(gs[1, 3:6])
axU = [fig.add_subplot(gs[2, 2 * i:2 * i + 2]) for i in range(3)]


def dot(ax, grps, xlabels, ds_names=True, seps=True, xpad=0.6):
    ncol = len(DS_ORDER) * len(grps)
    for gi, gname in enumerate(GENES):
        for di, key in enumerate(DS_ORDER):
            for pi, grp in enumerate(grps):
                rr = dp[(dp.dataset == key) & (dp.gene == gname) & (dp.group == grp)].iloc[0]
                xpos = di * len(grps) + pi
                if np.isnan(rr["mean"]):
                    if pi == max(0, len(grps) // 2 - (len(grps) == 2)):
                        ax.text(xpos + (0.5 if len(grps) == 2 else 0), gi, "n.m.",
                                ha="center", va="center", fontsize=7, color="#b8b6b0")
                    continue
                ax.scatter(xpos, gi, s=8 + 130 * rr.frac, c=[seq_blue(rr.scaled)],
                           edgecolors="#52514e", linewidths=0.4)
    ax.set_xticks(range(ncol)); ax.set_xticklabels(xlabels, fontsize=7.0)
    if seps:
        for di in range(1, len(DS_ORDER)):
            ax.axvline(di * len(grps) - 0.5, color="#e1e0d9", lw=0.8)
    if ds_names:
        for di, key in enumerate(DS_ORDER):
            ax.text(di * len(grps) + (len(grps) - 1) / 2, -1.15, DSL[key].split(" (")[0],
                    ha="center", fontsize=8.5, fontweight="bold")
    ax.set_ylim(len(GENES) - 0.4, -1.75); ax.set_xlim(-0.6, ncol - 1 + xpad)
    ax.tick_params(length=0)


# ---------- a: signature defined without FC ----------
dot(axA, ["BG", "OtherCBastro"], ["BG", "parenchymal\nastrocytes"] * 3)
axA.set_yticks(range(len(GENES)))
axA.set_yticklabels([f"$\\it{{{g}}}$" for g in GENES], fontsize=8)
axA.set_title("Bergmann signature — defined without FC", fontsize=9, loc="left", pad=14)
panel_letter(axA, "a", dx=-0.14, dy=1.10)

# ---------- b: FC independently reproduce it ----------
dot(axB1, ["FC"], ["", "", ""], ds_names=False, seps=False, xpad=1.0)
axB1.set_yticks(range(len(GENES))); axB1.set_yticklabels([])
axB1.set_title("FC\n(MERFISH · Stereo-seq · Xenium)", fontsize=7.5, loc="left", pad=6)
panel_letter(axB1, "b", dx=-0.16, dy=1.10)
axB2.scatter(xv[ok], yv[ok], s=22, c=GROUP["FC"], edgecolors="#7a1213", linewidths=0.4, zorder=5)
lim = [min(xv[ok].min(), yv[ok].min()) - 0.3, max(xv[ok].max(), yv[ok].max()) + 0.3]
axB2.plot(lim, lim, c="#898781", ls="--", lw=0.9, zorder=1)
axB2.axhline(0, c="#e1e0d9", lw=0.75); axB2.axvline(0, c="#e1e0d9", lw=0.75)
axB2.set_xlim(lim); axB2.set_ylim(lim)
axB2.set_xlabel("Δlog2, BG vs parenchymal\n(signature definition — no FC)", fontsize=7.5)
axB2.set_ylabel("Δlog2, FC vs parenchymal", fontsize=7.2, labelpad=2)
axB2.text(0.04, 0.97, f"r = {r:.2f} ({int(ok.sum())} gene × dataset pairs)\ny: independent FC contrast",
          transform=axB2.transAxes, fontsize=7.2, va="top")

# ---------- c: quantitative identity (kNN + centroid ratios) ----------
knn = pd.read_csv(TAB + "/fig2_knn_fc_in_bergmann.csv").set_index("key").loc[DS_ORDER]
assert knn.n_FC.astype(int).tolist() == [362, 89, 44]   # 4 animals (ABCA-1 integrated 2026-08-18)
assert np.allclose(knn.frac_BG, [98.9871086, 68.0898876, 88.4848485], atol=1e-4)
assert np.allclose(knn.base_BG, [76.1070800, 61.3526570, 50.2480706], atol=1e-4)
x = np.arange(3); w = 0.36
axC1.bar(x - w / 2, knn.frac_BG, w, color=GROUP["BG"], label="observed: neighbours that are BG")
axC1.bar(x + w / 2, knn.base_BG, w, color="#b8b6b0", label="chance: BG share of astroglia")
axC1.set_xticks(x); axC1.set_xticklabels([DSL[k].replace(" (", "\n(") for k in DS_ORDER], fontsize=8)
axC1.set_ylabel("% of FC's 15 nearest\nastroglial neighbours")
axC1.set_ylim(0, 132)
axC1.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0.01, 1.00))
panel_letter(axC1, "c", dx=-0.20, dy=1.12)

boot = pd.read_csv(TAB + "/fig_bootstrap_ratios.csv").set_index("dataset").loc[DS_ORDER]
assert np.allclose(boot.pooled, [3.5551720, 1.1470085, 2.1055181], atol=1e-4)
assert np.allclose(boot.med, [3.1464870, 1.1318223, 1.8953412], atol=1e-4)
sec = pd.read_csv(TAB + "/centroid_ratio_by_section.csv")
sec = sec[sec.level == "section"]
zh, ha = sec[sec.dataset == "zhuang"], sec[sec.dataset == "hao"]
assert len(zh) == 29 and (zh.ratio > 1).all(), "MERFISH per-section replication broken"
assert len(ha) == 14, len(ha)
ypos = dict(zip(DS_ORDER, np.arange(3)[::-1]))
rng_j = np.random.default_rng(0)
for key, ss in [("zhuang", zh), ("hao", ha)]:
    jit = rng_j.uniform(-0.17, 0.17, len(ss))
    axC2.scatter(ss.ratio, ypos[key] + jit, s=16, facecolors="none",
                 edgecolors="#52514e", linewidths=0.8, zorder=4,
                 label="per-section ratio" if key == "zhuang" else None)
axC2.errorbar(boot.loc["suyama", "pooled"], ypos["suyama"],
              xerr=[[boot.loc["suyama", "pooled"] - boot.lo["suyama"]],
                    [boot.hi["suyama"] - boot.loc["suyama", "pooled"]]],
              fmt="none", ecolor="#898781", elinewidth=1.0, capsize=2.5, zorder=3,
              label="cell bootstrap 95% (descriptive)")
axC2.vlines([boot.loc[k, "pooled"] for k in DS_ORDER], [ypos[k] - 0.26 for k in DS_ORDER],
            [ypos[k] + 0.26 for k in DS_ORDER], color="#0b0b0b", lw=1.8, zorder=5,
            label="pooled, all cells")
axC2.axvline(1.0, color="#898781", ls="--", lw=0.9)
axC2.set_yticks(list(ypos.values()))
axC2.set_yticklabels(["MERFISH\n29 sections", "Stereo-seq\n14 sections", "Xenium\n1 section"], fontsize=7.5)
axC2.set_xlabel("centroid ratio d(FC, parench.) / d(FC, BG)\n(Euclidean, full z-space; >1 = closer to BG)")
axC2.legend(frameon=False, fontsize=7, loc="upper right", bbox_to_anchor=(1.01, 1.04),
            handlelength=1.2, labelspacing=0.45)
axC2.set_xlim(0.72, 3.0); axC2.set_ylim(-0.85, 4.1)

# ---------- d: UMAPs (display projection, fixed seed) ----------
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
    ax.set_xlabel("UMAP1" if HAVE_UMAP else "PC1", fontsize=7)
    ax.set_ylabel("UMAP2" if HAVE_UMAP else "PC2", fontsize=7)
    ax.set_title(DSL[key], fontsize=9, fontweight="bold", pad=3)
    if i == 0:
        panel_letter(ax, "d", dx=-0.12, dy=1.16)

save(fig, "Figure3_identity")
