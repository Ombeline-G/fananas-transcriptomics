#!/usr/bin/env python
"""Figure S5 — FC-only reclustering: stability measured, not declared (review point 29).

The former version selected k by KMeans silhouette with ad hoc thresholds (0.15 declare /
0.25 benchmark). Replaced by bootstrap-consensus stability with a size-matched positive
control, per Antonin's recommendation ("show that no stable cluster exists"):
  - per dataset, FC-only: KMeans k = 2-4 refit on 100 random 80% subsamples (varying seeds);
    stability = mean pairwise adjusted Rand index (ARI) + consensus matrices (PAC released);
  - positive control at the SAME n: a 50/50 mix of canonical BG and parenchymal astrocytes
    of the same dataset, subsampled to the FC count (10 draws x 50 runs), identical pipeline —
    this measures what a genuine class boundary yields at exactly these sample sizes,
    making the power limitation (44-362 FC) explicit rather than implicit.
Result (4 animals, 2026-08-18): Stereo-seq, no reproducible partition (ARI <= 0.11 vs control
0.89); MERFISH, unstable (0.30 vs 0.99) — the full-data k=2 minor grouping (41/321) tracks
capture + Purkinje contamination; Xenium, ONE stable bipartition (0.86) = the technical quality axis
(31/13; 11/13 minor-state cells contamination-flagged; median detected genes 11 vs 20).
Stability confirms a partition, not a biology: the analysis is powered enough to find the
quality axis, and finds nothing else. Reported as a coarse check only (underpowered for
subtle states by design).

Exports: FigureS5_stability.csv, FigureS5_membership.csv.
(The former FigureS5_subcluster_DE / _section_composition / _subcluster_summary outputs
are retired with the threshold machinery.)
"""
import json, os, sys
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, SEMANTIC, GENE_PAL, DSL, DS_ORDER, W2, panel_letter, save
use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT); TAB = BASE + "/comparison/tables"
B_FC, B_CTRL, N_DRAWS, FRAC, KS = 100, 50, 10, 0.8, (2, 3, 4)
DS_COL = dict(zip(DS_ORDER, GENE_PAL[:3]))


def pcs_of(E):
    Z = ((E - E.mean()) / (E.std() + 1e-9)).fillna(0).values
    return PCA(min(20, Z.shape[1] - 1, Z.shape[0] - 1), random_state=0).fit_transform(Z)


def stability(X, k, B, seed0=0):
    """Mean pairwise ARI over B KMeans runs on 80% subsamples, + consensus matrix and PAC."""
    n = len(X); rng = np.random.default_rng(seed0)
    labs, idxs = [], []
    for b in range(B):
        idx = rng.choice(n, max(k * 3, int(FRAC * n)), replace=False)
        lab = KMeans(k, n_init=10, random_state=seed0 * 1000 + b).fit_predict(X[idx])
        labs.append(dict(zip(idx, lab))); idxs.append(set(idx))
    aris = []
    for i in range(B):
        for j in range(i + 1, min(i + 26, B)):
            common = list(idxs[i] & idxs[j])
            aris.append(adjusted_rand_score([labs[i][c] for c in common],
                                            [labs[j][c] for c in common]))
    co = np.zeros((n, n)); tot = np.zeros((n, n))
    for b in range(B):
        idx = np.array(sorted(idxs[b])); lab = np.array([labs[b][c] for c in idx])
        co[np.ix_(idx, idx)] += (lab[:, None] == lab[None, :]).astype(float)
        tot[np.ix_(idx, idx)] += 1
    cons = np.divide(co, tot, out=np.full_like(co, np.nan), where=tot > 0)
    iu = np.triu_indices(n, 1); vals = cons[iu][np.isfinite(cons[iu])]
    pac = float(((vals > 0.1) & (vals < 0.9)).mean())
    return float(np.mean(aris)), pac, cons


stab_rows, meta = [], {}
for key in DS_ORDER:
    m = pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    if "cell" not in m.columns: m = m.rename(columns={"index": "cell"})
    E = pd.read_parquet(f"{TAB}/{key}_astro_expr.parquet")
    m = m.set_index("cell"); E = E.reindex(m.index)
    fcm = m[m.group == "FC"]; Efc = E.reindex(fcm.index); nfc = len(fcm)
    X = pcs_of(Efc)
    cons2 = None
    for k in KS:
        a, p, cons = stability(X, k, B_FC, seed0=0)
        if k == 2: cons2 = cons
        stab_rows.append(dict(dataset=DSL[key], cohort="FC", k=k, n_cells=nfc,
                              ARI=round(a, 3), PAC=round(p, 3), n_draws=1))
        print(f"{key:7s} FC       n={nfc:3d} k={k}  ARI={a:.3f}  PAC={p:.3f}")
    # size-matched positive control: 50/50 BG + parenchymal, identical pipeline
    bg_i = m.index[m.group == "BG"]; pa_i = m.index[m.group == "OtherCBastro"]
    rng = np.random.default_rng(1); acc = {k: [] for k in KS}
    for d in range(N_DRAWS):
        nb = nfc // 2; na = nfc - nb
        sel = list(rng.choice(bg_i, nb, replace=False)) + list(rng.choice(pa_i, na, replace=False))
        Xc = pcs_of(E.reindex(sel))
        for k in KS:
            a, p, _ = stability(Xc, k, B_CTRL, seed0=d + 2)
            acc[k].append((a, p))
    for k in KS:
        arr = np.array(acc[k], float)
        stab_rows.append(dict(dataset=DSL[key], cohort="control_BG+parenchymal", k=k, n_cells=nfc,
                              ARI=round(arr[:, 0].mean(), 3), PAC=round(arr[:, 1].mean(), 3),
                              n_draws=N_DRAWS))
        print(f"{key:7s} control  n={nfc:3d} k={k}  ARI={arr[:,0].mean():.3f}  PAC={arr[:,1].mean():.3f}")
    lab2 = KMeans(2, n_init=10, random_state=0).fit_predict(X)   # full-data bipartition (descriptive)
    meta[key] = dict(fcm=fcm, Efc=Efc, X=X, lab2=lab2, cons2=cons2, nfc=nfc)

STAB = pd.DataFrame(stab_rows)
STAB.to_csv(f"{TAB}/FigureS5_stability.csv", index=False)

# ---- membership export + Xenium/Zhuang full-data k=2 characterization ----
flags = set(map(str, json.load(open(
    f"{PROJECT_ROOT}/datasets/suyama_2026_xenium/results/tables/suyama_fc_contaminated.json"))))
mem_rows = []
for key in DS_ORDER:
    D = meta[key]; nd = (D["Efc"] > 0).sum(1).values
    for c, l, n in zip(D["fcm"].index, D["lab2"], nd):
        mem_rows.append(dict(dataset=DSL[key], cell=c, k2_state=f"state-{l+1}", n_detected=int(n),
                             contamination_flag=bool(str(c) in flags) if key == "suyama" else False))
pd.DataFrame(mem_rows).to_csv(f"{TAB}/FigureS5_membership.csv", index=False)

# ================= AUDIT (G1) =================
S = STAB.set_index(["dataset", "cohort", "k"])
CANON = {("zhuang", "FC", 2): 0.302, ("hao", "FC", 2): 0.086, ("suyama", "FC", 2): 0.857,   # 4 animals (2026-08-18)
         ("zhuang", "control_BG+parenchymal", 2): 0.991,
         ("hao", "control_BG+parenchymal", 2): 0.886,
         ("suyama", "control_BG+parenchymal", 2): 0.623}
for (key, coh, k), v in CANON.items():
    got = float(S.loc[(DSL[key], coh, k), "ARI"])
    assert abs(got - v) < 5e-4, (key, coh, k, got)
    print(f"AUDIT ARI {key:7s} {coh:24s} k={k}: {got:.3f} (canon {v:.3f}) ✓")
assert float(S.loc[(DSL["hao"], "FC", 3), "ARI"]) < 0.15 and float(S.loc[(DSL["hao"], "FC", 4), "ARI"]) < 0.15
zh = meta["zhuang"]; sizes_zh = np.bincount(zh["lab2"])
assert sorted(sizes_zh) == [41, 321], sizes_zh                      # 4 animals: unstable minor grouping
nd_zh = (zh["Efc"] > 0).sum(1).values; mm_zh = zh["lab2"] == np.argmin(sizes_zh)
assert (np.median(nd_zh[mm_zh]), np.median(nd_zh[~mm_zh])) == (116.0, 90.0)
p_pk = mannwhitneyu(zh["fcm"]["contam_Purkinje"].values[mm_zh],
                    zh["fcm"]["contam_Purkinje"].values[~mm_zh])[1]
assert p_pk < 1e-4    # the minor grouping tracks capture + Purkinje contamination (and is unstable: ARI 0.30)
su = meta["suyama"]; sizes_su = np.bincount(su["lab2"])
assert sorted(sizes_su) == [13, 31], sizes_su
minor = int(np.argmin(sizes_su)); mm = su["lab2"] == minor
nflag = int(sum(str(c) in flags for c in su["fcm"].index[mm]))
nd_su = (su["Efc"] > 0).sum(1).values
assert nflag == 11, nflag                                            # 11/13 minor-state cells flagged
md_min, md_maj = float(np.median(nd_su[mm])), float(np.median(nd_su[~mm]))
assert (md_min, md_maj) == (11.0, 20.0), (md_min, md_maj)
p_nd = mannwhitneyu(nd_su[mm], nd_su[~mm])[1]
assert p_nd < 1e-4
print(f"AUDIT zhuang k=2 = 321/41, minor = capture+Purkinje axis (116 vs 90 genes, p_pk={p_pk:.1e}), unstable ✓")
print(f"AUDIT suyama k=2 = 31/13, minor: {nflag}/13 flagged, detected genes {md_min:.0f} vs {md_maj:.0f} (p={p_nd:.1e}) ✓")

# ================= FIGURE =================
fig = plt.figure(figsize=(W2, 5.6))
gs = GridSpec(2, 3, figure=fig, height_ratios=[1.0, 0.92], hspace=0.62, wspace=0.42,
              top=0.92, bottom=0.09, left=0.09, right=0.97)

# (a) stability curves: FC vs size-matched known-class control
axa = fig.add_subplot(gs[0, 0:2])
for key in DS_ORDER:
    fc = STAB[(STAB.dataset == DSL[key]) & (STAB.cohort == "FC")]
    ct = STAB[(STAB.dataset == DSL[key]) & (STAB.cohort != "FC")]
    axa.plot(fc.k, fc.ARI, "-o", color=DS_COL[key], ms=4, lw=1.2)
    axa.plot(ct.k, ct.ARI, "--s", color=DS_COL[key], ms=4, lw=1.0, mfc="white")
axa.set_xticks(list(KS)); axa.set_ylim(0, 1.02); axa.set_yticks([0, 0.5, 1])
axa.set_xlabel("k (KMeans)", fontsize=8); axa.set_ylabel("bootstrap stability (mean ARI)", fontsize=8)
axa.set_title("FC-only clustering never reaches its size-matched\nknown-class control (except the Xenium quality axis)",
              loc="left", fontsize=9, fontweight="bold")
handles = ([Line2D([], [], color=DS_COL[k], lw=1.2, label=DSL[k].split(" ")[0]) for k in DS_ORDER] +
           [Line2D([], [], color="#52514e", lw=1.2, marker="o", ms=4, label="FC only"),
            Line2D([], [], color="#52514e", lw=1.0, ls="--", marker="s", ms=4, mfc="white",
                   label="BG+parenchymal mix (same n)")])
axa.legend(handles=handles, fontsize=7, frameon=False, loc="upper center",
           bbox_to_anchor=(0.5, -0.30), ncol=3, columnspacing=1.1, handletextpad=0.5)
panel_letter(axa, "a")

# (c) the one stable partition (Xenium): the technical quality axis
axc = fig.add_subplot(gs[0, 2])
Xs = su["X"]; C_STATE = {0: GROUP["FC"], 1: SEMANTIC["contam"]}
order_states = [int(np.argmax(sizes_su)), minor]
for rank, st in enumerate(order_states):
    mm2 = su["lab2"] == st
    axc.scatter(Xs[mm2, 0], Xs[mm2, 1], s=16, c=C_STATE[rank], linewidths=0,
                label=f"state {rank+1} (n={int(mm2.sum())})")
isflag = np.array([str(c) in flags for c in su["fcm"].index])
axc.scatter(Xs[isflag, 0], Xs[isflag, 1], s=26, marker="x", c="#0b0b0b", linewidths=0.9,
            label="contamination flag")
axc.set_xticks([]); axc.set_yticks([]); axc.margins(0.10)
axc.set_xlabel("PC1", fontsize=8); axc.set_ylabel("PC2", fontsize=8)
axc.set_title("Xenium: the stable split is\nthe quality axis (11/13 flagged;\n11 vs 20 genes detected)",
              loc="left", fontsize=9, fontweight="bold")
axc.legend(fontsize=7, frameon=False, loc="upper left", handletextpad=0.2, labelspacing=0.35,
           borderaxespad=0.1)
panel_letter(axc, "c")

# (b) consensus matrices, k = 2, FC only — "no stable cluster" made visible
for i, key in enumerate(DS_ORDER):
    ax = fig.add_subplot(gs[1, i])
    cons = meta[key]["cons2"].copy(); np.fill_diagonal(cons, 1.0)
    dist = 1 - np.nan_to_num(cons, nan=0.5)
    order = leaves_list(linkage(dist[np.triu_indices(len(dist), 1)], method="average"))
    ax.imshow(cons[np.ix_(order, order)], cmap="Greys", vmin=0, vmax=1, aspect="equal",
              interpolation="nearest")
    a2 = float(S.loc[(DSL[key], "FC", 2), "ARI"])
    ax.set_title(f"{DSL[key].split(' ')[0]} — ARI {a2:.2f}", loc="left", fontsize=8.5)
    ax.set_xticks([]); ax.set_yticks([])
    if i == 0:
        ax.set_ylabel("FC cells (consensus-ordered)", fontsize=7.2)
        panel_letter(ax, "b")
fig.text(0.09, 0.015, "b: co-clustering frequency over 100 subsampled KMeans runs (k = 2; white 0 — black 1). "
         "MERFISH: an unstable minor grouping (41/321, ARI 0.30 vs control 0.99) tracking capture and Purkinje contamination;\n"
         "Stereo-seq: diffuse (no stable partition at any k); Xenium: two blocks — the quality axis shown in c.",
         fontsize=7, color="#52514e", va="bottom")
save(fig, "FigureS5_heterogeneity")
plt.close(fig)
