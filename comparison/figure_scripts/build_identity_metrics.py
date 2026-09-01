#!/usr/bin/env python3
"""Rebuild the identity metrics used by Figures 3, 4 and S3.

All expression matrices are z-scored within dataset. kNN uses 30 PCs and asks,
for each FC, how many of its 15 nearest non-FC astroglial neighbours are BG.

Centroid analysis (fully specified — review point 6):
- space: the full per-gene z-scored panel (all panel genes, no PCA);
- distance: Euclidean;
- centroids: per-group mean over ALL cells of that group (all BG, all
  parenchymal astrocytes, all FC of the dataset);
- ratio: d(FC, parenchymal) / d(FC, BG), >1 = FC closer to BG.
The cell-level bootstrap (1,000 resamples of cells with replacement within
each of the three groups jointly, percentile 2.5-97.5, independent fixed seed
per dataset) is DESCRIPTIVE within-dataset resampling uncertainty, not
biological uncertainty: cells of one animal/section are correlated. Resampling
noise inflates both centroid distances and therefore attenuates the ratio
toward 1, so the bootstrap summary is conservative for the claim ratio > 1
(its median sits below the plug-in pooled ratio, column `pooled`).
Replication is instead assessed section by section (centroid_ratio_by_section
.csv): within each section holding >=MIN_FC_SECTION FC, >=10 BG and >=10
parenchymal astrocytes, all three centroids are recomputed from that section's
cells only (same global z-space) — this removes section-level batch shifts and
treats sections, and animals for multi-animal datasets, as replication units.
Few-FC centroid noise attenuates these ratios toward 1 as well (conservative);
`ratio_noise_corrected` subtracts the expected sampling-noise term tr(S)/n
from each squared distance and confirms per-section ratios consistent with the
pooled value. Undefined for the single-section Xenium dataset.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.neighbors import NearestNeighbors

from project_paths import ROOT

TABLES = ROOT / "comparison" / "tables"
ORDER = ["zhuang", "hao", "suyama"]
NAMES = {
    "zhuang": "MERFISH",
    "hao": "Stereo-seq",
    "suyama": "Xenium",
}
YEAR_NAMES = {
    "zhuang": "MERFISH",
    "hao": "Stereo-seq",
    "suyama": "Xenium",
}
SEEDS = {"zhuang": 0, "hao": 1, "suyama": 2}
GROUPS = ["FC", "BG", "OtherCBastro"]
MIN_FC_SECTION = 3
MIN_REF_SECTION = 10


def centroid_ratio(z, groups, sel):
    """Plug-in ratio d(FC, parenchymal)/d(FC, BG) on cells `sel`, plus the
    noise-corrected variant (squared distances minus tr(S)/n per group;
    NaN when a corrected squared distance is not positive)."""
    stats = {}
    for group in GROUPS:
        X = z[sel & (groups == group)]
        stats[group] = (X.mean(axis=0), X.var(axis=0, ddof=1).sum() / len(X))
    d_fc_bg = float(np.linalg.norm(stats["FC"][0] - stats["BG"][0]))
    d_fc_other = float(np.linalg.norm(stats["FC"][0] - stats["OtherCBastro"][0]))

    def corrected_sq(a, b):
        return ((stats[a][0] - stats[b][0]) ** 2).sum() - stats[a][1] - stats[b][1]

    c_bg, c_other = corrected_sq("FC", "BG"), corrected_sq("FC", "OtherCBastro")
    corrected = float(np.sqrt(c_other / c_bg)) if c_bg > 0 and c_other > 0 else float("nan")
    return d_fc_bg, d_fc_other, d_fc_other / d_fc_bg, corrected


def load(key):
    meta = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    if "cell" not in meta.columns:
        meta = meta.rename(columns={"index": "cell"})
    expr = pd.read_parquet(TABLES / f"{key}_astro_expr.parquet")
    if meta.cell.isin(expr.index).all():
        expr = expr.reindex(meta.cell)
    assert len(meta) == len(expr)
    z = ((expr - expr.mean()) / (expr.std() + 1e-9)).fillna(0).to_numpy()
    return meta, z


knn_rows, distance_rows, bootstrap_rows, section_rows = [], [], [], []
silhouette_cells = []
for key in ORDER:
    meta, z = load(key)
    groups = meta.group.to_numpy()
    fc = np.flatnonzero(groups == "FC")
    reference = np.flatnonzero(groups != "FC")
    assert len(fc) and len(reference) >= 15

    n_components = min(30, z.shape[0] - 1, z.shape[1])
    pcs = PCA(n_components=n_components, random_state=0).fit_transform(z)
    neighbours = NearestNeighbors(n_neighbors=15).fit(pcs[reference]).kneighbors(
        pcs[fc], return_distance=False
    )
    neighbour_groups = groups[reference[neighbours]]
    frac_bg = 100 * float((neighbour_groups == "BG").mean())
    n_bg = int((groups == "BG").sum())
    n_other = int((groups == "OtherCBastro").sum())
    base_bg = 100 * n_bg / (n_bg + n_other)
    knn_rows.append(dict(name=NAMES[key], key=key, frac_BG=frac_bg, base_BG=base_bg, n_FC=len(fc)))

    centroids = {group: z[groups == group].mean(axis=0) for group in GROUPS}
    d_fc_other = float(np.linalg.norm(centroids["FC"] - centroids["OtherCBastro"]))
    d_bg_fc = float(np.linalg.norm(centroids["BG"] - centroids["FC"]))
    d_bg_other = float(np.linalg.norm(centroids["BG"] - centroids["OtherCBastro"]))
    # Pairwise two-group silhouettes (review point 8): FC|BG against the two
    # accepted-class references. Computed in the same 30-PC space as the kNN
    # analysis, where class separation is visible; the full z-scored panel is
    # kept for the record but distance concentration in >=1,000 noisy
    # dimensions compresses ALL silhouettes toward 0 there (even BG vs
    # parenchymal), which made the old 3-group panel unreadable.
    silhouette_rows = [("silhouette_3grp", float(silhouette_score(z, groups)))]
    for a, b in [("FC", "BG"), ("BG", "OtherCBastro"), ("FC", "OtherCBastro")]:
        sel = (groups == a) | (groups == b)
        samples = silhouette_samples(pcs[sel], groups[sel])
        silhouette_rows += [
            (f"silhouette2_{a}_{b}_pc30", float(samples.mean())),
            (f"silhouette2_{a}_{b}_full", float(silhouette_score(z[sel], groups[sel]))),
        ]
        silhouette_cells.append(pd.DataFrame({
            "dataset": key, "pair": f"{a}|{b}", "group": groups[sel],
            "s_pc30": samples.astype("float32")}))
    for pair, value in [
        ("FC-OtherCBastro", d_fc_other),
        ("BG-FC", d_bg_fc),
        ("BG-OtherCBastro", d_bg_other),
    ] + silhouette_rows:
        distance_rows.append(dict(dataset=YEAR_NAMES[key], pair=pair, dist=value))
    pooled_ratio = d_fc_other / d_bg_fc

    rng = np.random.default_rng(SEEDS[key])
    indices = {group: np.flatnonzero(groups == group) for group in GROUPS}
    ratios = []
    for _ in range(1000):
        boot_centroids = {
            group: z[rng.choice(idx, len(idx), replace=True)].mean(axis=0)
            for group, idx in indices.items()
        }
        d_fc_bg = np.linalg.norm(boot_centroids["FC"] - boot_centroids["BG"])
        d_fc_other = np.linalg.norm(boot_centroids["FC"] - boot_centroids["OtherCBastro"])
        ratios.append(d_fc_other / d_fc_bg)
    bootstrap_rows.append(dict(
        dataset=key,
        pooled=float(pooled_ratio),
        med=float(np.median(ratios)),
        lo=float(np.percentile(ratios, 2.5)),
        hi=float(np.percentile(ratios, 97.5)),
    ))

    # --- per-section / per-animal replication of the centroid ratio (review point 6) ---
    animals = meta.animal.to_numpy()
    sections = meta.section.to_numpy()
    units = [("section", s, (sections == s)) for s in pd.unique(sections)]
    if meta.animal.nunique() > 1:
        units += [("animal", a, (animals == a)) for a in pd.unique(animals)]
    if len(pd.unique(sections)) > 1:
        for level, unit, sel in units:
            n = {group: int((groups[sel] == group).sum()) for group in GROUPS}
            if n["FC"] < MIN_FC_SECTION or n["BG"] < MIN_REF_SECTION or n["OtherCBastro"] < MIN_REF_SECTION:
                continue
            d_bg, d_other, ratio, corrected = centroid_ratio(z, groups, sel)
            section_rows.append(dict(
                dataset=key, level=level, unit=unit,
                animal=animals[sel][0] if level == "section" else unit,
                n_FC=n["FC"], n_BG=n["BG"], n_parenchymal=n["OtherCBastro"],
                d_FC_BG=d_bg, d_FC_parenchymal=d_other,
                ratio=ratio, ratio_noise_corrected=corrected,
            ))

pd.DataFrame(knn_rows).to_csv(TABLES / "fig2_knn_fc_in_bergmann.csv", index=False)
pd.DataFrame(distance_rows).to_csv(TABLES / "fig4_distance.csv", index=False)
pd.DataFrame(bootstrap_rows).to_csv(TABLES / "fig_bootstrap_ratios.csv", index=False)
sec = pd.DataFrame(section_rows)
sec.to_csv(TABLES / "centroid_ratio_by_section.csv", index=False)
pd.concat(silhouette_cells).to_csv(TABLES / "silhouette_cells_pc30.csv", index=False)
print(pd.DataFrame(knn_rows).to_string(index=False))
print(pd.DataFrame(bootstrap_rows).to_string(index=False))
for (key, level), grp in sec.groupby(["dataset", "level"]):
    above = int((grp.ratio > 1).sum())
    print(f"{key} per-{level}: n={len(grp)}  median {grp.ratio.median():.2f}  "
          f"range {grp.ratio.min():.2f}-{grp.ratio.max():.2f}  >1: {above}/{len(grp)}  "
          f"noise-corrected median {grp.ratio_noise_corrected.median():.2f}")
