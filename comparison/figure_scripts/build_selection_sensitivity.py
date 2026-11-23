#!/usr/bin/env python3
"""Sensitivity of the main results to the FC selection margin and to manual review.

Reviewer request (S. Alon, November 2026): (i) do the conclusions on the molecular
composition of FC hold at other margins than 45 µm above the PCL line (e.g. 35, 55 µm)?
(ii) what do the results look like WITHOUT the manual rejection step?

MERFISH (four animals) — the selection margin is the segmentation parameter
ML_MIN_DEPTH (layer_segmentation/segment_cerebellar_layers.py). The per-cell
assignments store every Bergmann-class soma of the ML mask with its distance to the
PCL, so the geometric candidate set can be rebuilt exactly for any margin D without
re-running the segmentation: ML Bergmann somata with distance >= D, then the same
declumping rule (> DECLUMP_NB other candidates within DECLUMP_R -> removed). At
D = 45 µm this reproduces the swipe candidate pool exactly (asserted: 240 + 406).

Arms
- validated : manually validated FC whose gating distance is >= D (D >= 35). The cohort of
              the study is D = 45 (swipe decisions fc_decisions_zhuang.csv + fc_decisions_abca1.csv);
              the 35-µm cohort adds the candidates of the 35-45 µm band accepted in the dedicated
              review (swipe_band35.html -> swipe_validation/fc_decisions_band35.csv, same rater,
              same criteria, blind to expression).
- geometric : every geometric candidate at margin D, NO manual review (D >= 35).
- rejected  : the candidates of the 45-µm pool discarded by manual review, on their own.

BG (confident PCL Bergmann) and parenchymal astrocytes are the canonical sets and
are identical in every arm; only the FC set changes. Every readout re-uses the
canonical definition (build_identity_metrics, build_equivalence, harmonized_lean,
build_full_de); the validated arm at 45 µm must reproduce the canonical numbers.

Xenium — the anti-PCL rule is a per-soma distance to the nearest Purkinje soma
(canonical >= 30 µm); stricter variants are subsets of the canonical 44 FC.
Stereo-seq — validated bins by depth (10-µm grid; canonical minimum 50 µm).
For these two platforms only the validated arm is available from the active tables.

Output: tables/selection_sensitivity.csv (one row per dataset x arm x margin) and
tables/selection_sensitivity_DE_hits.csv (every capture-matched DE hit per MERFISH arm).
Cache: tables/zhuang_sensitivity_extra_expr.parquet (expression of the geometric
candidates that are not in the canonical cohort; rebuilt if missing).
"""
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests

from panels import module_scores
from project_paths import ROOT, TABLES, ZHUANG, SUYAMA

MARGINS_GEOM = [35, 45, 55, 65, 80, 100]
MARGINS_VALID = [45, 55, 65, 80, 100]
DECLUMP_R, DECLUMP_NB = 0.035, 2          # mm; identical to segment_cerebellar_layers.py
GROUPS = ["FC", "BG", "OtherCBastro"]
MODULE = {"zhuang": ["Aqp4", "Agt", "Gja1", "Aldh1l1"],
          "hao": ["Agt", "Slc1a2", "Slc1a3", "S100b", "Gjb6"],
          "suyama": ["Aqp4", "Agt", "Gja1", "Slc1a3", "S100b", "Aldh1l1", "Gjb6"]}
CALIPER = 0.10


# ------------------------------------------------------------------ candidates
def ml_bergmann(seg):
    """Every Bergmann-class soma of the ML mask, before any margin or declumping."""
    amb = seg.region == "Ambiguous_BG"
    return seg[(seg.region == "FC")
               | (amb & (seg.dist_to_PCL_um >= 45))                              # declumped at 45
               | (amb & (seg.confidence > 0.30) & (seg.dist_to_PCL_um < 45))]    # ML, shallower


def geometric_candidates(seg, margin):
    m = ml_bergmann(seg)
    m = m[m.dist_to_PCL_um >= margin]
    keep = []
    for _, g in m.groupby("section"):
        xy = g[["x", "y"]].to_numpy()
        tree = cKDTree(xy)
        nb = np.array([len(tree.query_ball_point(p, DECLUMP_R)) - 1 for p in xy])
        keep += list(g.cell_label[nb <= DECLUMP_NB])
    return set(keep)


# ------------------------------------------------------------------ readouts
def match_pairs(meta, capture):
    pure = meta[(meta.contam_Purkinje < 1) & (meta.contam_granule < 1)]
    fc, bg = pure[pure.group == "FC"], pure[pure.group == "BG"]
    pairs = []
    for sec, f in fc.groupby("section"):
        b = bg[bg.section == sec]
        if b.empty:
            continue
        cap_b = capture[b.index]
        used = set()
        for i in f.index:
            d = (cap_b - capture[i]).abs()
            d = d[[j for j in d.index if j not in used]]
            if len(d) == 0:
                continue
            j = d.idxmin()
            if d[j] <= CALIPER * max(capture[i], 1):
                used.add(j)
                pairs.append((i, j))
    return pairs


def readouts(meta, expr, key, de=True, tost=True):
    """meta: index = cell, columns group/section/animal(+contam); expr aligned."""
    out = {}
    groups = meta.group.to_numpy()
    n = {g: int((groups == g).sum()) for g in GROUPS}
    out.update(n_FC=n["FC"], n_BG=n["BG"], n_parenchymal=n["OtherCBastro"])
    # --- identity: kNN, centroid ratio, silhouettes (build_identity_metrics) ---
    z = ((expr - expr.mean()) / (expr.std() + 1e-9)).fillna(0).to_numpy()
    fc, ref = np.flatnonzero(groups == "FC"), np.flatnonzero(groups != "FC")
    pcs = PCA(n_components=min(30, z.shape[0] - 1, z.shape[1]), random_state=0).fit_transform(z)
    nb = NearestNeighbors(n_neighbors=15).fit(pcs[ref]).kneighbors(pcs[fc], return_distance=False)
    out["knn_pct_BG"] = 100 * float((groups[ref[nb]] == "BG").mean())
    out["knn_base_pct_BG"] = 100 * n["BG"] / (n["BG"] + n["OtherCBastro"])
    cen = {g: z[groups == g].mean(axis=0) for g in GROUPS}
    out["centroid_ratio"] = float(np.linalg.norm(cen["FC"] - cen["OtherCBastro"])
                                  / np.linalg.norm(cen["FC"] - cen["BG"]))
    for a, b in [("FC", "BG"), ("BG", "OtherCBastro")]:
        sel = (groups == a) | (groups == b)
        out[f"silhouette_{a}_{b}".replace("OtherCBastro", "par")] = float(
            silhouette_samples(pcs[sel], groups[sel]).mean())
    # --- effect-size magnitudes (build_equivalence) ---
    det = pd.DataFrame({g: (expr[groups == g] > 0).mean() for g in GROUPS})
    E = expr.loc[:, det.max(axis=1) >= 0.10]
    mean = {g: E[groups == g].mean() for g in GROUPS}
    d_fc_bg, d_bg_par, d_fc_par = (mean["FC"] - mean["BG"], mean["BG"] - mean["OtherCBastro"],
                                   mean["FC"] - mean["OtherCBastro"])
    out.update(n_expressed_genes=E.shape[1],
               med_abs_d_FC_BG=float(d_fc_bg.abs().median()),
               med_abs_d_BG_par=float(d_bg_par.abs().median()),
               pct_genes_gt025_FC_BG=100 * float((d_fc_bg.abs() > 0.25).mean()),
               pct_genes_gt025_BG_par=100 * float((d_bg_par.abs() > 0.25).mean()),
               r_FCpar_BGpar=float(np.corrcoef(d_fc_par, d_bg_par)[0, 1]))
    if tost:
        counts = meta.groupby("section").group.value_counts().unstack(fill_value=0)
        counts = counts.reindex(columns=GROUPS, fill_value=0)
        secs = counts.index[(counts.FC >= 3) & (counts.BG >= 10) & (counts.OtherCBastro >= 10)]
        out["n_tost_sections"] = len(secs)
        if len(secs) >= 3:
            sec_arr = meta.section.to_numpy()
            per = pd.DataFrame([E[(sec_arr == s) & (groups == "FC")].mean()
                                - E[(sec_arr == s) & (groups == "BG")].mean() for s in secs])
            mu, se = per.mean(), per.std(ddof=1) / np.sqrt(len(secs))
            tc = stats.t.ppf(0.95, len(secs) - 1)
            for margin in (0.25, 0.5):
                eq = ((mu - tc * se) > -margin) & ((mu + tc * se) < margin)
                out[f"pct_equiv_{margin}"] = 100 * float(eq.mean())
    # --- astrocytic shift (harmonized_lean) ---
    zz = (expr - expr.mean()) / expr.std(ddof=0)
    capture = (expr > 0).sum(axis=1).astype(float)
    pairs = match_pairs(meta, capture)
    fi, bi = [p[0] for p in pairs], [p[1] for p in pairs]
    mod = [g for g in MODULE[key] if g in expr.columns]
    score = zz[mod].mean(axis=1)
    diffs = score[fi].to_numpy() - score[bi].to_numpy()
    out.update(n_pairs=len(pairs), lean_delta=float(np.mean(diffs)) if len(diffs) else np.nan,
               lean_p=float(wilcoxon(diffs).pvalue) if len(diffs) >= 6 else np.nan)
    if "animal" in meta.columns and meta.animal.nunique() > 1:
        per_animal = pd.Series(diffs, index=meta.loc[fi, "animal"].to_numpy()).groupby(level=0).mean()
        out["lean_by_animal"] = ";".join(f"{a.replace('Zhuang-', '')}:{v:+.3f}" for a, v in per_animal.items())
        out["lean_animals_positive"] = f"{int((per_animal > 0).sum())}/{len(per_animal)}"
    for g in ["Aqp4", "Gfap"]:
        if g in expr.columns and len(pairs) >= 6:
            gd = zz[g][fi].to_numpy() - zz[g][bi].to_numpy()
            out[f"dz_{g}"] = float(gd.mean())
            out[f"p_{g}"] = float(wilcoxon(gd).pvalue) if np.abs(gd).sum() > 0 else np.nan
    # --- capture-matched DE (build_full_de.paired_de) ---
    hits = pd.DataFrame()
    if de and len(pairs) >= 6:
        X, Y = expr.loc[fi].to_numpy(), expr.loc[bi].to_numpy()
        D = X - Y
        pv = np.ones(D.shape[1])
        for k in range(D.shape[1]):
            if np.abs(D[:, k]).sum() > 0:
                try:
                    pv[k] = wilcoxon(D[:, k]).pvalue
                except Exception:
                    pv[k] = 1.0
        padj = multipletests(pv, method="fdr_bh")[1]
        res = pd.DataFrame({"gene": expr.columns, "log2FC": D.mean(axis=0), "padj": padj})
        hits = res[(res.padj < 0.05) & (res.log2FC.abs() > 0.25)].sort_values("log2FC")
        out["n_DE_hits"] = len(hits)
        out["DE_hits"] = ";".join(f"{g}({v:+.2f})" for g, v in zip(hits.gene, hits.log2FC))
    return out, hits


# ================================================================== MERFISH
cells = pd.read_parquet(TABLES / "zhuang_astro_cells.parquet").set_index("cell")
expr = pd.read_parquet(TABLES / "zhuang_astro_expr.parquet").reindex(cells.index)
kept = pd.read_csv(ZHUANG / "layer_segmentation" / "BG_assignments_kept.csv", dtype={"cell_label": str})
a1 = pd.read_csv(ZHUANG / "layer_segmentation" / "BG_assignments_ABCA1.csv", dtype={"cell_label": str})
a1_swipe = pd.read_csv(ZHUANG / "layer_segmentation" / "FC_candidates_ABCA1_for_swipe.csv",
                       dtype={"cell_label": str})
a1 = a1[a1.section.isin(a1_swipe.section.unique())]          # the 14 swiped ABCA-1 sections
seg = pd.concat([kept, a1], ignore_index=True)
seg_d = seg.drop_duplicates("cell_label").set_index("cell_label")

dec = pd.concat([pd.read_csv(ZHUANG / "swipe_validation" / f, dtype={"cell_label": str})
                 for f in ["fc_decisions_zhuang.csv", "fc_decisions_abca1.csv"]])
reviewed = set(dec.cell_label)
validated = set(dec.cell_label[dec.decision == "yes"])
canon_fc = set(cells.index[cells.group == "FC"])
assert validated == canon_fc and len(validated) == 362
geom = {D: geometric_candidates(seg, D) for D in MARGINS_GEOM}
assert geom[45] == reviewed and len(reviewed) == 646, "45-µm reconstruction != swipe pool"
band_35_45 = {c for c in geom[35] if seg_d.dist_to_PCL_um[c] < 45}
print("geometric candidates per margin:", {D: len(s) for D, s in geom.items()},
      "| 35-45 band:", len(band_35_45))

# expression of candidates outside the canonical cohort
extra_ids = sorted(set().union(*geom.values()) - set(cells.index))
cache = TABLES / "zhuang_sensitivity_extra_expr.parquet"
if cache.exists() and set(pd.read_parquet(cache).index) >= set(extra_ids):
    extra = pd.read_parquet(cache)
else:
    import anndata as ad
    import scipy.sparse as sp
    frames, want = [], set(extra_ids)
    for an in ["Zhuang-ABCA-1", "Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4"]:
        f = ZHUANG / "abc_cache" / "expression_matrices" / an / "20230830" / f"{an}-log2.h5ad"
        a = ad.read_h5ad(f, backed="r")
        m = a.obs.index.astype(str).isin(want)
        if m.sum() == 0:
            continue
        sub = a[m].to_memory()
        X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
        frames.append(pd.DataFrame(X, columns=sub.var["gene_symbol"].astype(str).values,
                                   index=sub.obs.index.astype(str)))
        print("  loaded", an, int(m.sum()), flush=True)
    extra = pd.concat(frames)
    extra = extra.loc[~extra.index.duplicated()]
    extra.to_parquet(cache)
extra = extra.reindex(columns=expr.columns)
assert set(extra_ids) <= set(extra.index), "expression missing for some geometric candidates"
assert not (set(extra_ids) & set(cells.index))

ref = cells[cells.group != "FC"][["group", "section", "animal"]]
all_expr = pd.concat([expr, extra.loc[extra_ids]])


def merfish_arm(fc_ids):
    # canonical cells keep the row order of the canonical table (the greedy capture
    # matching is order-dependent), unreviewed/rejected candidates follow, sorted by id
    fc_ids = [c for c in cells.index if c in fc_ids] + sorted(set(fc_ids) - set(cells.index))
    sec = seg_d.section.reindex(fc_ids)
    fc = pd.DataFrame({"group": "FC", "section": sec.values,
                       "animal": [s.rsplit(".", 1)[0] for s in sec.values]}, index=fc_ids)
    meta = pd.concat([fc, ref])
    E = all_expr.loc[meta.index]
    sco, _ = module_scores(E)
    meta[["contam_Purkinje", "contam_granule"]] = sco[["contam_Purkinje", "contam_granule"]].to_numpy()
    return meta, E


rows, hit_rows = [], []


def run(dataset, arm, margin, meta, E, key, **kw):
    out, hits = readouts(meta, E, key, **kw)
    row = dict(dataset=dataset, arm=arm, margin_um=margin, **out)
    row["FC_share_pct"] = 100 * out["n_FC"] / (out["n_FC"] + out["n_BG"])
    rows.append(row)
    for g, v, p in zip(hits.get("gene", []), hits.get("log2FC", []), hits.get("padj", [])):
        hit_rows.append(dict(dataset=dataset, arm=arm, margin_um=margin, gene=g, log2FC=v, padj=p))
    print(f"{dataset:10s} {arm:15s} >={margin:>3} µm  n_FC={out['n_FC']:4d}  kNN={out['knn_pct_BG']:.1f}%  "
          f"ratio={out['centroid_ratio']:.2f}  sil={out['silhouette_FC_BG']:+.3f}  "
          f"med|d|={out['med_abs_d_FC_BG']:.3f}  lean={out['lean_delta']:+.3f} (p={out['lean_p']:.1e}, "
          f"{out['n_pairs']} pairs)  DE={out.get('n_DE_hits', 'NA')}", flush=True)


gate_d = seg_d.dist_to_PCL_um
# 35-45 µm band: reviewed separately (swipe_band35.html -> fc_decisions_band35.csv); the accepted
# candidates extend the validated arm to 35 µm. The cohort of the study (>= 45 µm) is unchanged.
band_cand = ZHUANG / "layer_segmentation" / "FC_candidates_band35_45_for_swipe.csv"
assert set(pd.read_csv(band_cand, dtype={"cell_label": str}).cell_label) == band_35_45, \
    "swipe candidate file of the 35-45 µm band != band rebuilt here"
bd = pd.read_csv(ZHUANG / "swipe_validation" / "fc_decisions_band35.csv", dtype={"cell_label": str})
assert set(bd.cell_label) == band_35_45 and not bd.cell_label.duplicated().any(), \
    "band decisions do not cover exactly the 35-45 µm candidates"
assert bd.decision.isin(["yes", "no"]).all(), "unresolved decisions (maybe/empty) in fc_decisions_band35.csv"
band_yes = set(bd.cell_label[bd.decision == "yes"])
print(f"35-45 µm band reviewed: {len(band_yes)}/{len(bd)} accepted")
run("MERFISH", "validated", 35, *merfish_arm(validated | band_yes), "zhuang")
for D in MARGINS_VALID:
    run("MERFISH", "validated", D, *merfish_arm({c for c in validated if gate_d[c] >= D}), "zhuang")
for D in MARGINS_GEOM:
    run("MERFISH", "geometric", D, *merfish_arm(geom[D]), "zhuang")
# rejected-only cohort: the candidates the reviewer discarded, as a direct look at what review removes
run("MERFISH", "rejected", 45, *merfish_arm(reviewed - validated), "zhuang")

# ================================================================== Xenium / Stereo-seq (validated subsets)
def simple(key):
    c = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    if "cell" in c.columns:
        c = c.set_index("cell")
    elif "index" in c.columns:
        c = c.set_index("index")
    c.index = c.index.astype(str)
    e = pd.read_parquet(TABLES / f"{key}_astro_expr.parquet")
    e.index = c.index
    return c, e


sc, se = simple("suyama")
strict = set(map(str, json.load(open(SUYAMA / "results" / "tables" / "suyama_fc_gate_strict45.json"))))
assert strict == set(sc.index[(sc.group == "FC") & (sc.depth_um >= 45)]), "depth_um != d_pkj of the strict set"
hc, he = simple("hao")
# margin = depth above the PCL (Stereo-seq, 10-µm grid) or distance to the nearest Purkinje soma (Xenium)
for name, key, c, e, margins in [("Stereo-seq", "hao", hc, he, [45, 55, 65, 80, 100]),
                                 ("Xenium", "suyama", sc, se, [30, 35, 40, 45])]:
    for D in margins:
        keep = (c.group != "FC") | (c.depth_um >= D)
        if (keep & (c.group == "FC")).sum() < 10:
            continue
        run(name, "validated", D, c[keep], e[keep.to_numpy()], key, de=False, tost=False)

out = pd.DataFrame(rows)
out.to_csv(TABLES / "selection_sensitivity.csv", index=False)
pd.DataFrame(hit_rows).to_csv(TABLES / "selection_sensitivity_DE_hits.csv", index=False)
print("saved selection_sensitivity.csv", out.shape)
