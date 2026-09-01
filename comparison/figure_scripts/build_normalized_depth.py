#!/usr/bin/env python3
"""Normalized ML depth of every FC (review point 4, Antonin): 0 = PCL, 1 = pia.

Per-dataset estimator (documented in Methods):
  Zhuang  — exact per-soma ratio depth/(depth + dist_to_pia) from the census geometry
            (FC matched to zhuang_ML_census_cells by cell id).
  Hao     — depth / local ML thickness. Local thickness = 95th percentile of the
            PCL-distance of ML-annotated bins within 100 um of the FC; per-bin PCL
            distance = distance to the nearest segmented BG soma (same PCL proxy as
            the pipeline's dist_to_PCL_um). Coordinates x10 um.
  Suyama  — |sdist| / local ML thickness, sdist = signed distance to the granule-mass
            boundary (suyama_geom_layers.csv); local thickness = 95th percentile of
            |sdist| of ML-band cells within 100 um of the FC.
Output: tables/fc_normalized_depth.csv  (dataset, cell, norm_depth, method)
"""
import os, sys
import numpy as np, pandas as pd, anndata as ad
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"
HAO = BASE + "/datasets/hao_2024_stereo_seq"
SU = BASE + "/datasets/suyama_2026_xenium"
R_LOCAL = 100.0   # um, lateral neighbourhood for the local-thickness estimate
PCT = 95

out = []

# ---------------- Zhuang: exact ----------------
fc = pd.read_parquet(TAB + "/zhuang_astro_cells.parquet")
fc = fc[fc.group == "FC"].set_index("cell")
cen = pd.read_parquet(TAB + "/zhuang_ML_census_cells.parquet").set_index("cell")
common = fc.index.intersection(cen.index)
nd = cen.loc[common, "depth_from_PCL"] / (cen.loc[common, "depth_from_PCL"] + cen.loc[common, "dist_to_pia"])
for c, v in nd.items():
    out.append(dict(dataset="zhuang", cell=c, norm_depth=float(v), method="exact depth/(depth+pia)"))
print(f"Zhuang: {len(nd)}/{len(fc)} FC matched to census geometry | median {nd.median():.3f} (IQR {nd.quantile(.25):.3f}-{nd.quantile(.75):.3f})")

# ---------------- Zhuang ABCA-1: local ML thickness from the segmentation masks ----------------
# The census (and its pia geometry) covers ABCA-2/3/4 only; ABCA-1 FC (integrated 2026-08-18)
# use the same local-thickness estimator as Hao/Suyama, computed from the section masks:
# thickness = PCT-th percentile of the PCL-distance of ML-mask pixels within R_LOCAL of the FC.
import json
ZH = BASE + "/datasets/zhang_zhuang_2023_merfish"
a1 = fc[fc.animal == "Zhuang-ABCA-1"]
nd1 = []
for sec, sub in a1.groupby("section"):
    sd = f"{ZH}/layer_segmentation/sections/{sec}"
    if not os.path.exists(sd + "/grid_meta.json"):
        continue
    meta = json.load(open(sd + "/grid_meta.json"))
    ml_m = np.load(sd + "/ML_mask.npy"); pcl_m = np.load(sd + "/PCL_mask.npy")
    um = meta["res"] * 1000.0
    ml_px = np.argwhere(ml_m)[:, ::-1] * um            # (x, y) in um on the grid frame
    pcl_px = np.argwhere(pcl_m)[:, ::-1] * um
    if len(ml_px) < 20 or len(pcl_px) < 5:
        continue
    d_ml_pcl, _ = cKDTree(pcl_px).query(ml_px)
    tree_ml = cKDTree(ml_px)
    for c, row in sub.iterrows():
        px = ((row.x - meta["x0"]) * 1000.0, (row.y - meta["y0"]) * 1000.0)
        idx = tree_ml.query_ball_point(px, r=R_LOCAL)
        if len(idx) < 10 or not np.isfinite(row.depth_um):
            continue
        thick = np.percentile(d_ml_pcl[idx], PCT)
        if thick <= 0:
            continue
        v = float(row.depth_um / thick)
        nd1.append(v)
        out.append(dict(dataset="zhuang", cell=c, norm_depth=v, method="depth / local ML thickness (masks)"))
nd1 = pd.Series(nd1)
print(f"Zhuang ABCA-1: {len(nd1)}/{len(a1)} FC | median {nd1.median():.3f} (IQR {nd1.quantile(.25):.3f}-{nd1.quantile(.75):.3f})")

# ---------------- Hao: local ML thickness from native annotation ----------------
seg = pd.read_csv(HAO + "/layer_segmentation/BG_assignments_all.csv")
seg = seg[seg.section.astype(str).str.startswith("Mouse1_")]
fcs = seg[seg.region == "FC"].copy()
dec = pd.read_csv(HAO + "/swipe_validation/fc_decisions_hao.csv")
dec = dec[(dec.animal == "Mouse1") & (dec.decision == "yes")]
import re
strip = lambda s: re.sub(r"-\d+$", "", str(s))
valid = set(zip(dec.section, dec.cell_label.map(strip)))
fcs["cl"] = fcs.cell_label.map(strip)
fcs = fcs[[ (s, c) in valid for s, c in zip(fcs.section, fcs.cl) ]]
UM = 10.0
vals = []
for sec, sub in fcs.groupby("section"):
    p = f"{HAO}/data_raw/stereo_seq_mouse1/{sec}.h5ad"
    if not os.path.exists(p):
        continue
    a = ad.read_h5ad(p, backed="r")
    ob = a.obs
    if "annotation" not in ob.columns:
        continue
    ml = ob[ob.annotation.astype(str) == "molecular layer"]
    if not hasattr(ml, "x") and "x" not in ml.columns:
        # coordinates from obsm if not in obs
        xy = pd.DataFrame(a.obsm["spatial"], index=ob.index, columns=["x", "y"]).loc[ml.index]
    else:
        xy = ml[["x", "y"]]
    bgsec = seg[(seg.section == sec) & (seg.region == "BG")]
    if len(bgsec) < 5 or len(xy) < 20:
        continue
    tree_bg = cKDTree(bgsec[["x", "y"]].values * UM)
    mlc = xy.values * UM
    dpcl_ml, _ = tree_bg.query(mlc)
    tree_ml = cKDTree(mlc)
    for _, rrow in sub.iterrows():
        pos = np.array([rrow.x, rrow.y]) * UM
        idx = tree_ml.query_ball_point(pos, R_LOCAL)
        if len(idx) < 20:
            continue
        T = np.percentile(dpcl_ml[idx], PCT)
        if T <= 0:
            continue
        v = min(rrow.dist_to_PCL_um / T, 1.05)
        vals.append(v)
        out.append(dict(dataset="hao", cell=f"{sec}|{rrow.cl}", norm_depth=float(v),
                        method=f"depth / local ML thickness (p{PCT}, r={R_LOCAL:.0f}um)"))
vals = np.array(vals)
print(f"Hao: {len(vals)}/{len(fcs)} FC with local thickness | median {np.median(vals):.3f} (IQR {np.percentile(vals,25):.3f}-{np.percentile(vals,75):.3f})")

# ---------------- Suyama: local ML thickness from geometric band ----------------
gl = pd.read_csv(SU + "/results/tables/suyama_geom_layers.csv")
gl["cell"] = gl["cell"].astype(str)
cells = pd.read_parquet(TAB + "/suyama_astro_cells.parquet")
sfc = cells[cells.group == "FC"].copy(); sfc["cell"] = sfc["cell"].astype(str)
g = gl.set_index("cell")
mlg = gl[gl.layer_geom == "ML"]
tree_ml = cKDTree(mlg[["x", "y"]].values)
sd_ml = np.abs(mlg.sdist_um.values)
vals = []
for _, rrow in sfc.iterrows():
    if rrow.cell not in g.index:
        continue
    srow = g.loc[rrow.cell]
    pos = np.array([srow.x, srow.y]) if np.isfinite(srow.x) else None
    if pos is None:
        continue
    idx = tree_ml.query_ball_point(pos, R_LOCAL)
    if len(idx) < 20:
        continue
    T = np.percentile(sd_ml[idx], PCT)
    if T <= 0:
        continue
    v = min(abs(srow.sdist_um) / T, 1.05)
    vals.append(v)
    out.append(dict(dataset="suyama", cell=rrow.cell, norm_depth=float(v),
                    method=f"|sdist| / local ML thickness (p{PCT}, r={R_LOCAL:.0f}um)"))
vals = np.array(vals)
print(f"Suyama: {len(vals)}/{len(sfc)} FC | median {np.median(vals):.3f} (IQR {np.percentile(vals,25):.3f}-{np.percentile(vals,75):.3f})")

pd.DataFrame(out).to_csv(TAB + "/fc_normalized_depth.csv", index=False)
print("saved fc_normalized_depth.csv")
