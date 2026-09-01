#!/usr/bin/env python
"""Figure S4 — the cross-platform signature table (ONE heatmap-table; review point 28).

Rows: genes significant (padj < 0.05, |Δlog2| > 0.25) with concordant direction in >= 2
datasets in the BG-vs-parenchymal contrast — the contrast that DEFINES the signature
(no FC cell enters row admission; contamination-panel genes excluded).
Columns: 3 contrasts x 3 datasets (DS_ORDER) —
  (1) BG vs parenchymal astrocytes   — defines the signature;
  (2) FC vs parenchymal astrocytes   — the independent test;
  (3) FC vs BG, capture-matched      — READ from the canonical Figure-4 DE tables
      (fig4_DE_*_depthmatched.csv), never recomputed here. The former in-script
      FC-vs-BG contrast matched FC and BG on soma depth (a spatial coordinate) and
      contradicted the Methods rule that position is never matched; it is removed.

Exports (stable names):
  fig6_all_contrasts_DE.csv    — BG_vs_OtherAstro + FC_vs_OtherAstro, all genes
                                 (consumed by Figure3_identity.py, r = 0.95 scatter)
  Figure6_signature_table.csv  — wide per-gene table: the figure's exact data
                                 (FigureS3_confounds.py appends purkinje_spatial_covar)
"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panels import PANELS
from figstyle import use_style, DSL, DS_ORDER, W2, save
use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT); TAB = BASE + "/comparison/tables"
CONTAM_GENES = {g for cat, gs in PANELS.items() if cat.startswith("contam_") for g in gs}
SIG_TH, FDR = 0.25, 0.05


def de(A, B):
    r = []
    for g in A.columns:
        x, y = A[g].values, B[g].values
        if len(x) < 3 or len(y) < 3:
            p = 1.0; d = np.nan
        else:
            try: _, p = mannwhitneyu(x, y, alternative="two-sided")
            except Exception: p = 1.0
            d = float(x.mean() - y.mean())
        r.append(dict(gene=g, log2FC=d, pval=p))
    out = pd.DataFrame(r)
    out["padj"] = multipletests(out.pval.fillna(1), method="fdr_bh")[1]
    return out


# ---- contrasts vs parenchymal astrocytes (within dataset; BG defines, FC tests) ----
allde = []
for key in DS_ORDER:
    m = pd.read_parquet(f"{TAB}/{key}_astro_cells.parquet")
    if "cell" not in m.columns: m = m.rename(columns={"index": "cell"})
    E = pd.read_parquet(f"{TAB}/{key}_astro_expr.parquet")
    m = m.set_index("cell"); E = E.reindex(m.index); grp = m.group
    for cname, gA in (("BG_vs_OtherAstro", "BG"), ("FC_vs_OtherAstro", "FC")):
        d = de(E[(grp == gA).values], E[(grp == "OtherCBastro").values])
        d["dataset"] = DSL[key]; d["contrast"] = cname; allde.append(d)
DE = pd.concat(allde, ignore_index=True)
DE.to_csv(f"{TAB}/fig6_all_contrasts_DE.csv", index=False)

# ---- FC vs BG: canonical capture-matched tables of the Figure-4 analysis ----
FCBG = {key: pd.read_csv(f"{TAB}/fig4_DE_{key}_depthmatched.csv").set_index("gene")
        for key in DS_ORDER}

# ---- row admission: BG-vs-parenchymal contrast only (non-circular) ----
sub = DE[DE.contrast == "BG_vs_OtherAstro"]
rows = []
for g, gg in sub.groupby("gene"):
    gg = gg.dropna(subset=["log2FC"])
    up = gg[(gg.padj < FDR) & (gg.log2FC > SIG_TH)]
    dn = gg[(gg.padj < FDR) & (gg.log2FC < -SIG_TH)]
    n_sig = max(len(up), len(dn))
    rows.append(dict(gene=g, n_present=gg.dataset.nunique(), n_sig=n_sig,
                     direction="Bergmann class" if len(up) >= len(dn) else "parenchymal astrocytes",
                     mean_BG=gg.log2FC.mean()))
AGG = pd.DataFrame(rows)
KEEP = (AGG[(AGG.n_sig >= 2) & (~AGG.gene.isin(CONTAM_GENES))]
        .sort_values("mean_BG", ascending=False).reset_index(drop=True))
GENES = KEEP.gene.tolist()
N_BERG = int((KEEP.direction == "Bergmann class").sum())


def fcbg_cell(g, key):
    """(Δlog2, significant?) for one gene in one dataset's canonical FC-vs-BG table."""
    t = FCBG[key]
    if g not in t.index: return np.nan, False
    v, q = float(t.loc[g, "log2FC"]), float(t.loc[g, "padj"])
    return v, (q < FDR and abs(v) > SIG_TH)


def block(contrast):
    """genes x datasets matrices of Δlog2 and significance for a parenchymal contrast."""
    s = DE[DE.contrast == contrast].set_index(["gene", "dataset"])
    V = np.full((len(GENES), 3), np.nan); S = np.zeros((len(GENES), 3), bool)
    for i, g in enumerate(GENES):
        for j, key in enumerate(DS_ORDER):
            ix = (g, DSL[key])
            if ix in s.index:
                V[i, j] = s.loc[ix, "log2FC"]
                S[i, j] = (s.loc[ix, "padj"] < FDR) and (abs(V[i, j]) > SIG_TH)
    return V, S


V_BG, S_BG = block("BG_vs_OtherAstro")
V_FC, S_FC = block("FC_vs_OtherAstro")
V_D = np.full((len(GENES), 3), np.nan); S_D = np.zeros((len(GENES), 3), bool)
for i, g in enumerate(GENES):
    for j, key in enumerate(DS_ORDER):
        V_D[i, j], S_D[i, j] = fcbg_cell(g, key)

# ================= AUDIT (G1) — every displayed value vs canonical =================
# constants recomputed 2026-08-18 after ABCA-1 integration (4 MERFISH animals)
CANON_BG = {"Hopx": (2.3404, 3), "Gpr37l1": (1.4349, 3), "A2m": (1.2724, 3), "Agt": (-2.1382, 3),
            "Gdf10": (2.7546, 2), "Slc1a3": (2.5135, 2), "Timp4": (1.7305, 2),
            "Fabp7": (1.6904, 2), "Tnc": (1.1410, 2)}
CANON_FC = {"Hopx": 2.2967, "Gpr37l1": 1.0886, "Agt": -2.0830,
            "A2m": 1.1739, "Gdf10": 3.0790, "Slc1a3": 2.1854, "Timp4": 1.5683}
for g, (mv, ns) in CANON_BG.items():
    got = float(AGG.set_index("gene").loc[g, "mean_BG"]); nsg = int(AGG.set_index("gene").loc[g, "n_sig"])
    assert abs(got - mv) < 5e-4 and nsg == ns, (g, got, nsg)
    print(f"AUDIT BG  {g:8s} mean {got:+.4f} (canon {mv:+.4f}) n_sig {nsg} ✓")
fcm = DE[DE.contrast == "FC_vs_OtherAstro"].groupby("gene").log2FC.mean()
for g, mv in CANON_FC.items():
    assert abs(float(fcm[g]) - mv) < 5e-4, (g, fcm[g])
    print(f"AUDIT FC  {g:8s} mean {float(fcm[g]):+.4f} (canon {mv:+.4f}) ✓")
assert len(GENES) == 19 and N_BERG == 12, (len(GENES), N_BERG)
assert not set(GENES) & CONTAM_GENES
zh_sig = {GENES[i] for i in range(len(GENES)) if S_D[i, 0]}
assert zh_sig == {"Tnc", "Ntsr2", "Gfap", "Aqp4", "Agt"}, zh_sig       # the 5 signature genes among the 13 MERFISH hits
assert not any(S_D[:, 1]) and not any(S_D[:, 2])                        # Hao (0 canonical), Xenium (0 hits)
assert [int(np.isfinite(V_BG[:, j]).sum()) for j in range(3)] == [13, 17, 19]   # MERFISH panel lacks 6; Hao curated panel lacks Homer3, Ednrb
assert [int(np.isfinite(V_D[:, j]).sum()) for j in range(3)] == [13, 12, 19]    # Hao pipeline gate excludes 7
print(f"AUDIT rows {len(GENES)} (12 Bergmann-side / 7 parenchymal-side), "
      f"FC-vs-BG sig: MERFISH {sorted(zh_sig)}, Stereo-seq none, Xenium none ✓")

# ---- wide signature table (the figure's exact data; FigureS3 appends covar flag) ----
SHORT = {"zhuang": "MERFISH", "hao": "StereoSeq", "suyama": "Xenium"}
wide = dict(Gene=GENES,
            Enriched_in=KEEP.direction.tolist(),
            n_datasets_sig=KEEP.n_sig.astype(int).tolist())
for tag, (V, S) in (("BGvsAstro", (V_BG, S_BG)), ("FCvsAstro", (V_FC, S_FC)), ("FCvsBG", (V_D, S_D))):
    for j, key in enumerate(DS_ORDER):
        wide[f"{tag}_{SHORT[key]}"] = np.round(V[:, j], 3)
    wide[f"sig_{tag}"] = [";".join(SHORT[DS_ORDER[j]] for j in range(3) if S[i, j]) or ""
                          for i in range(len(GENES))]
pd.DataFrame(wide).to_csv(f"{TAB}/Figure6_signature_table.csv", index=False)
print("signature table rows:", len(GENES))

# ================= FIGURE — one heatmap-table =================
VMAX = 2.5
BLOCKS = [("BG vs parenchymal\n(defines the signature)", V_BG, S_BG),
          ("FC vs parenchymal\n(independent test)", V_FC, S_FC),
          ("FC vs BG\n(capture-matched)", V_D, S_D)]
fig, axes = plt.subplots(1, 3, figsize=(W2, 7.1), gridspec_kw=dict(wspace=0.10))
fig.subplots_adjust(left=0.155, right=0.865, top=0.865, bottom=0.055)
CMAP = plt.cm.RdBu_r.copy(); CMAP.set_bad("#eceef1")   # n.m. cells light grey, distinct from white = 0
for b, (ax, (title, V, S)) in enumerate(zip(axes, BLOCKS)):
    im = ax.imshow(np.ma.masked_invalid(np.clip(V, -VMAX, VMAX)), cmap=CMAP, vmin=-VMAX, vmax=VMAX, aspect="auto")
    for i in range(V.shape[0]):
        for j in range(3):
            if np.isnan(V[i, j]):
                ax.text(j, i, "n.m.", ha="center", va="center", fontsize=7, color="#898781")
            else:
                ax.text(j, i, f"{V[i, j]:+.2f}{'*' if S[i, j] else ''}", ha="center",
                        va="center", fontsize=7,
                        color="#0b0b0b" if abs(V[i, j]) < 0.60 * VMAX else "white")
    ax.axhline(N_BERG - 0.5, color="#0b0b0b", lw=1.0)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["MERFISH", "Stereo-seq", "Xenium"], fontsize=8, rotation=45, ha="right")
    ax.set_title(title, fontsize=8.5, pad=5)
    if b == 0:
        ax.set_yticks(range(len(GENES)))
        ax.set_yticklabels([f"$\\it{{{g}}}$" for g in GENES], fontsize=8)
    else:
        ax.set_yticks([])
axes[0].text(-0.62, (N_BERG - 1) / 2, "higher in Bergmann class",
             transform=axes[0].get_yaxis_transform(), rotation=90, va="center",
             fontsize=8, fontweight="bold", color="#52514e")
axes[0].text(-0.62, N_BERG + (len(GENES) - N_BERG - 1) / 2, "higher in\nparenchymal",
             transform=axes[0].get_yaxis_transform(), rotation=90, va="center",
             fontsize=8, fontweight="bold", color="#52514e", ha="center")
cax = fig.add_axes([0.885, 0.35, 0.016, 0.30])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Δlog2 (clipped at ±2.5)", fontsize=8)
fig.text(0.155, 0.965, "Cross-platform signature table", fontsize=9.5, fontweight="bold")
fig.text(0.155, 0.935, "rows: significant, same direction in ≥2 datasets, BG-vs-parenchymal contrast only\n"
         "* padj < 0.05 and |Δlog2| > 0.25 within dataset · n.m. = not measured / not tested",
         fontsize=7.2, color="#52514e", va="top")
save(fig, "FigureS4_signatures")
plt.close(fig)
