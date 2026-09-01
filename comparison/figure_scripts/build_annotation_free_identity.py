#!/usr/bin/env python3
"""Annotation-independent identity test of ML astroglia (review 3f, design Antonin, 2026-08-17).

Training uses ONLY anatomically unambiguous reference cells (no ML cell):
  BG-PCL (canonical BG cohort), Astro-NT inside the WM, Astro-CB inside the GCL
  (robust_layers_zhuang.parquet x Allen subclass), n=800 each, seed 0.
Test set: the 519 ML-census astroglia, labels withheld.
Eight classifications: {multinomial logistic regression, random forest(400),
kNN-15, nearest correlation-centroid} x {full 1,122-gene panel, contamination-free panel}.
Signatures: BG-PCL vs Astro-NT-WM (50 genes each way, contamination-free space, no ML cells).

Outputs (comparison/tables/):
  annotfree_scores_cells.csv    per-cell bg/nt scores + majority prediction + positions
  annotfree_summary.csv         % assigned per census group x run
  annotfree_sig_genes.csv       the two signature gene lists
"""
import sys, os
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"
ZH = BASE + "/datasets/zhang_zhuang_2023_merfish"
rng = np.random.default_rng(0)
N = 800
CLASSES = ["BG", "Astro-NT", "Astro-CB"]
CONTAM = set(("Calb1 Car8 Pcp2 Pcp4 Itpr1 Grid2 Aldoc Slc1a6 Homer3 Gabra6 Gng13 Neurod1 "
              "Cbln1 Cbln3 Snap25 Stmn2 Syn1 Sncb Slc17a7 Fat2 Rbfox3 Adcy1 Nrxn3 "
              "Mbp Plp1 Mog Mag Cnp Mobp Cldn5 Pecam1 Rgs5 Pdgfrb Flt1 Dcn Col1a1 Igf2 Vtn Malat1").split())

# ---- reference and test cells ----
r = pd.read_parquet(TAB + "/robust_layers_zhuang.parquet")
mv = pd.read_parquet(ZH + "/cb_cells_multivote.parquet")
j = r.join(mv[["subclass"]], how="left"); scb = j.subclass.astype(str)
ntwm = j[(scb.str.contains("Astro-NT")) & (j.layer_robust == "WM")].index.astype(str)
cbgcl = j[(scb.str.contains("Astro-CB")) & (j.layer_robust == "GCL")].index.astype(str)
refc = pd.read_parquet(TAB + "/zhuang_astro_cells.parquet").set_index("cell")
bg = refc.index[refc.group == "BG"].astype(str)
cen = pd.read_parquet(TAB + "/zhuang_ML_census_cells.parquet")
cen = cen[cen.depth_from_PCL >= 45].copy()
cen["dist_border"] = cen[["dist_to_WM", "dist_to_pia"]].min(axis=1)
cen["grp4"] = np.where(cen.type == "Bergmann", "ML-Bergmann",
              np.where(cen.type == "Astro-CB", "ML-CB",
              np.where(cen.dist_border <= 30, "ML-NT-border", "ML-NT-interior")))
bg_tr = pd.Index(rng.choice(bg, N, replace=False))
nt_tr = pd.Index(rng.choice(ntwm, N, replace=False))
cb_tr = pd.Index(rng.choice(cbgcl, N, replace=False))
need = set(nt_tr) | set(cb_tr) | set(cen.cell.astype(str))
frames = []
for an in ["Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4"]:
    a = ad.read_h5ad(f"{ZH}/abc_cache/expression_matrices/{an}/20230830/{an}-log2.h5ad", backed="r")
    m = a.obs.index.astype(str).isin(need)
    if not m.sum():
        continue
    sub = a[m].to_memory()
    X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    frames.append(pd.DataFrame(X, columns=sub.var["gene_symbol"].astype(str).values,
                               index=sub.obs.index.astype(str)))
FE = pd.concat(frames); FE = FE.loc[~FE.index.duplicated()]
re_ = pd.read_parquet(TAB + "/zhuang_astro_expr.parquet"); re_.index = refc.index.astype(str)
genes = [g for g in re_.columns if g in FE.columns]
E_bg = re_.loc[bg_tr, genes]
E_nt = FE.reindex(nt_tr)[genes].dropna()
E_cb = FE.reindex(cb_tr)[genes].dropna()
cen = cen[cen.cell.astype(str).isin(FE.index)]
E_ml = FE.reindex(cen.cell.astype(str))[genes]
print(f"training BG {len(E_bg)} / NT-WM {len(E_nt)} / CB-GCL {len(E_cb)} | test {len(E_ml)} | genes {len(genes)}")

# ---- eight classifications ----
pred_cols = {}
for tag, gset in [("full", genes), ("clean", [g for g in genes if g not in CONTAM])]:
    Xtr = pd.concat([E_bg[gset], E_nt[gset], E_cb[gset]])
    ytr = np.array(["BG"] * len(E_bg) + ["Astro-NT"] * len(E_nt) + ["Astro-CB"] * len(E_cb))
    mu, sd = Xtr.mean(), Xtr.std() + 1e-9
    Ztr = ((Xtr - mu) / sd).values
    Zml = ((E_ml[gset] - mu) / sd).values
    pred_cols[f"logreg_{tag}"] = LogisticRegression(max_iter=3000, class_weight="balanced").fit(Ztr, ytr).predict(Zml)
    pred_cols[f"rforest_{tag}"] = RandomForestClassifier(400, random_state=0, class_weight="balanced").fit(Ztr, ytr).predict(Zml)
    pred_cols[f"kNN15_{tag}"] = KNeighborsClassifier(15).fit(Ztr, ytr).predict(Zml)
    cent = {c: Ztr[ytr == c].mean(0) for c in CLASSES}
    cor = np.array([[np.corrcoef(z, cent[c])[0, 1] for c in CLASSES] for z in Zml])
    pred_cols[f"corrcent_{tag}"] = np.array(CLASSES)[cor.argmax(1)]
P = pd.DataFrame(pred_cols, index=E_ml.index)
maj = P.apply(lambda row: row.value_counts().idxmax(), axis=1)

# ---- summary table ----
rows = []
for grp in ["ML-Bergmann", "ML-NT-border", "ML-NT-interior", "ML-CB"]:
    m = (cen["grp4"] == grp).values
    for run in P.columns:
        for c in CLASSES:
            rows.append(dict(group=grp, n=int(m.sum()), run=run, assigned=c,
                             pct=round(100 * float((P.loc[m.values if False else P.index[m], run] == c).mean()), 1)))
summary = pd.DataFrame(rows)
summary.to_csv(TAB + "/annotfree_summary.csv", index=False)

# ---- signatures (contamination-free space, BG-PCL vs Astro-NT-WM only) ----
gset = [g for g in genes if g not in CONTAM]
tr = pd.concat([E_bg[gset], E_nt[gset]])
mu, sd = tr.mean(), tr.std()
keep = [g for g in gset if sd[g] > 1e-6]
mu, sd = mu[keep], sd[keep]
zb = (E_bg[keep] - mu) / sd
zn = (E_nt[keep] - mu) / sd
diff = zb.mean() - zn.mean()
pv = []
for g in keep:
    try:
        pv.append(mannwhitneyu(zb[g], zn[g])[1])
    except Exception:
        pv.append(1.0)
padj = multipletests(np.nan_to_num(np.array(pv), nan=1.0), method="fdr_bh")[1]
ok = pd.Series(padj < 0.05, index=keep)
BGsig = diff[ok].sort_values(ascending=False).head(50).index.tolist()
NTsig = diff[ok].sort_values().head(50).index.tolist()
pd.DataFrame({"signature": ["BG"] * 50 + ["AstroNT"] * 50, "gene": BGsig + NTsig}).to_csv(
    TAB + "/annotfree_sig_genes.csv", index=False)

def score(E, gl):
    return ((E[keep] - mu) / sd)[gl].mean(axis=1)

out = []
for name, E in [("ref BG-PCL", E_bg), ("ref AstroNT-WM", E_nt), ("ref CB-GCL", E_cb)]:
    out.append(pd.DataFrame({"set": name, "grp4": "", "bg_score": score(E, BGsig),
                             "nt_score": score(E, NTsig), "pred_majority": "", "x": np.nan,
                             "y": np.nan, "section": ""}))
cpos = cen.set_index(cen.cell.astype(str))
out.append(pd.DataFrame({"set": "ML-census", "grp4": cpos.grp4, "bg_score": score(E_ml, BGsig),
                         "nt_score": score(E_ml, NTsig), "pred_majority": maj,
                         "x": cpos.x, "y": cpos.y, "section": cpos.section}))
cells = pd.concat(out)
cells.index.name = "cell"
cells.to_csv(TAB + "/annotfree_scores_cells.csv")
print("majority predictions:", maj.groupby(cen["grp4"].values).value_counts().to_dict())
print("saved annotfree_{scores_cells,summary,sig_genes}.csv")
