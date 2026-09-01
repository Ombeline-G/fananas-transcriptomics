#!/usr/bin/env python3
"""Platform-scoping diagnostics for the equivalence analysis (review point 7).

Documents WHY the magnitude/equivalence analysis (build_equivalence.py) is
restricted to MERFISH, from the raw sources:

- Stereo-seq (Hao, Mouse1): the median |Delta-log2|(FC-BG) across the 2,573
  canonical astro bins is as large as or larger than the BG-parenchymal class
  contrast (ratio 1.39), persists after exact capture equalization
  (multinomial downsampling of every bin to 600 transcripts, the canonical
  depth-matched procedure, seed 0; ratio 1.42), and is carried by the Purkinje
  program (Car8, Calb1, Pcp2, Aldoc) — bin composition, not astroglial
  biology, dominates magnitudes.
- Xenium (Suyama): median |Delta-log2|(FC-BG) is at class-contrast level
  (ratio 0.93) and the largest magnitudes are neuronal and oligodendroglial
  transcripts (Snap25, Sncb, Camk2b, Cnp) — position-correlated neuropil
  spillover that a layer contrast cannot be separated from without the
  matched design already used for DE.
The corr(Delta, detection-difference) column is reported for completeness but
does NOT discriminate platforms (log-normalized Delta and detection are
mathematically coupled; MERFISH shows r = 0.84 with a median |Delta| of 0.03).
The discriminators are the magnitude of the bulk distribution, its
carrier genes, and its persistence under capture equalization.

Sources: Hao section h5ads (data_raw/stereo_seq_mouse1), Xenium MEX
(datasets/suyama_2026_xenium/data/cell_feature_matrix — extracted from the
Zenodo archive P0_P6_Adult_mousecerebellum_5K.tar.gz, whole-slide run).
Slow (~20 min): reads 32 h5ads and a 600 MB MEX matrix.

Output: tables/fig4_equivalence_platform_scoping.csv
"""
import gzip

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

from project_paths import ROOT, HAO, SUYAMA

TABLES = ROOT / "comparison" / "tables"
GROUPS = ["FC", "BG", "OtherCBastro"]
MIN_DETECTION = 0.10
THIN_TARGET = 600
rows = []


def contrasts(E, groups):
    """Median |pooled Delta-log2| per contrast on the expressed universe,
    plus corr(Delta FC-BG, detection difference) and top |Delta| genes."""
    detection = pd.DataFrame({g: (E[groups == g] > 0).mean() for g in GROUPS})
    Ee = E.loc[:, detection.max(axis=1) >= MIN_DETECTION]
    mean = {g: Ee[groups == g].mean() for g in GROUPS}
    d = pd.DataFrame({"FC_BG": mean["FC"] - mean["BG"],
                      "BG_par": mean["BG"] - mean["OtherCBastro"],
                      "FC_par": mean["FC"] - mean["OtherCBastro"]})
    det_diff = (detection.loc[Ee.columns, "FC"] - detection.loc[Ee.columns, "BG"])
    top = d.FC_BG.abs().sort_values(ascending=False).head(10)
    return dict(
        n_genes=Ee.shape[1],
        med_FC_BG=float(d.FC_BG.abs().median()),
        med_BG_par=float(d.BG_par.abs().median()),
        med_FC_par=float(d.FC_par.abs().median()),
        ratio=float(d.FC_BG.abs().median() / d.BG_par.abs().median()),
        corr_delta_detection=float(np.corrcoef(d.FC_BG, det_diff)[0, 1]),
        top_FC_BG_genes=";".join(top.index),
    )


# ---------------- Hao: raw and capture-equalized ----------------
meta = pd.read_parquet(TABLES / "hao_astro_cells.parquet").rename(columns={"index": "cell"})
targets = set(meta.cell)
frames = []
for sec in sorted(meta.section.unique()):
    a = ad.read_h5ad(HAO / "data_raw" / "stereo_seq_mouse1" / f"{sec}.h5ad")
    a.var_names_make_unique()
    cid = np.array([f"{sec}|{i}" for i in a.obs.index.astype(str)])
    keep = np.flatnonzero(np.isin(cid, list(targets)))
    X = a.X[keep]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    frames.append(pd.DataFrame(np.rint(X).astype(np.int64),
                               columns=list(map(str, a.var_names)), index=cid[keep]))
    print(sec, len(keep), flush=True)
C = pd.concat(frames).fillna(0).astype(np.int64)
C = C.loc[~C.index.duplicated()]
hao_meta = meta.set_index("cell").reindex(C.index)
tot = C.sum(axis=1).to_numpy().astype(float)
tot[tot == 0] = 1

raw = pd.DataFrame(np.log1p(C.to_numpy() / tot[:, None] * 1e4) / np.log(2),
                   columns=C.columns, index=C.index)
res = contrasts(raw, hao_meta.group.to_numpy())
rows.append(dict(platform="Stereo-seq", variant="raw bins", **res))

ok = tot >= THIN_TARGET
rng = np.random.default_rng(0)
M = C.to_numpy()
D = np.zeros((int(ok.sum()), M.shape[1]), dtype=np.int64)
for j, i in enumerate(np.flatnonzero(ok)):
    D[j] = rng.multinomial(THIN_TARGET, M[i] / tot[i])
thin = pd.DataFrame(np.log1p(D / THIN_TARGET * 1e4) / np.log(2),
                    columns=C.columns, index=C.index[ok])
res = contrasts(thin, hao_meta.group.to_numpy()[ok])
rows.append(dict(platform="Stereo-seq", variant=f"capture-equalized {THIN_TARGET}", **res))

# ---------------- Suyama: raw 5K panel ----------------
MEX = SUYAMA / "data" / "cell_feature_matrix"
feat = pd.read_csv(MEX / "features.tsv.gz", sep="\t", header=None, names=["id", "name", "type"])
with gzip.open(MEX / "barcodes.tsv.gz", "rt") as fh:
    barcodes = np.array([line.strip() for line in fh])
sm = pd.read_parquet(TABLES / "suyama_astro_cells.parquet")
if "cell" not in sm.columns:
    sm = sm.rename(columns={"index": "cell"})
sm["cell"] = sm.cell.astype(str)
hit = np.isin(barcodes, list(sm.cell))
assert hit.sum() == len(sm), (hit.sum(), len(sm))
print("reading Xenium MEX matrix ...", flush=True)
X = mmread(MEX / "matrix.mtx.gz").tocsc()
Xs = X[:, np.flatnonzero(hit)].toarray().T.astype(np.float64)
cells = barcodes[np.flatnonzero(hit)]
ge = (feat.type == "Gene Expression").to_numpy()
stot = Xs[:, ge].sum(1)
stot[stot == 0] = 1
E = pd.DataFrame(np.log1p(Xs[:, ge] / stot[:, None] * 1e4) / np.log(2),
                 columns=feat.name[ge].values, index=cells).reindex(sm.cell)
old = pd.read_parquet(TABLES / "suyama_astro_expr.parquet").reindex(sm.cell)
shared = [c for c in old.columns if c in E.columns]
assert np.abs(E[shared].to_numpy() - old[shared].to_numpy()).max() < 1e-3, \
    "MEX recomputation does not match the canonical 61-gene table"
res = contrasts(E, sm.group.to_numpy())
rows.append(dict(platform="Xenium", variant="raw cells", **res))

out = pd.DataFrame(rows)
out.to_csv(TABLES / "fig4_equivalence_platform_scoping.csv", index=False)
print(out.drop(columns="top_FC_BG_genes").to_string(index=False))
for _, r in out.iterrows():
    print(f"{r.platform} [{r.variant}] top |d| FC-BG: {r.top_FC_BG_genes}")
