#!/usr/bin/env python3
"""kNN sensitivity for the FC-in-Bergmann classification (added 2026-08-17, Antonin's request).

Replicates build_identity_metrics.py exactly (per-gene z over all astroglia ->
PCA(<=30 components, random_state=0) -> neighbours among non-FC astroglia), then varies:
  1. k sweep: k in {5, 10, 15, 30}, full neighbour pool.
  2. Leave-own-section-out (LOSO): for each FC the pool excludes its own section (k=15).
  3. Same-animal, cross-section: pool = same animal, other sections (k=15).
     Only informative for Zhuang (3 animals); FC whose restricted pool is smaller
     than k are skipped and counted.
Base rates are recomputed per variant as the mean BG share of each FC's allowed pool.
Output: tables/knn_sensitivity.csv
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from project_paths import TABLES

ORDER = ["zhuang", "hao", "suyama"]
K_SWEEP = [5, 10, 15, 30]
K_MAIN = 15


def load(key):
    meta = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    if "cell" not in meta.columns:
        meta = meta.rename(columns={"index": "cell"})
    expr = pd.read_parquet(TABLES / f"{key}_astro_expr.parquet")
    if meta.cell.isin(expr.index).all():
        expr = expr.reindex(meta.cell)
    assert len(meta) == len(expr)
    z = ((expr - expr.mean()) / (expr.std() + 1e-9)).fillna(0).to_numpy()
    return meta.reset_index(drop=True), z


def frac_bg_restricted(pcs, groups, fc_idx, pool_masks, k):
    """pool_masks: list of boolean masks over all cells, one per FC (its allowed pool)."""
    hits, tot, base_terms, used = 0, 0, [], 0
    for i, mask in zip(fc_idx, pool_masks):
        pool = np.flatnonzero(mask)
        if len(pool) < k:
            continue
        nn = NearestNeighbors(n_neighbors=k).fit(pcs[pool])
        nb = nn.kneighbors(pcs[[i]], return_distance=False)[0]
        lab = groups[pool[nb]]
        hits += int((lab == "BG").sum())
        tot += k
        base_terms.append(float((groups[pool] == "BG").mean()))
        used += 1
    if tot == 0:
        return np.nan, np.nan, 0, len(fc_idx)
    return 100 * hits / tot, 100 * float(np.mean(base_terms)), used, len(fc_idx) - used


rows = []
for key in ORDER:
    meta, z = load(key)
    groups = meta.group.to_numpy()
    sections = meta.section.astype(str).to_numpy()
    animals = meta.animal.astype(str).to_numpy() if "animal" in meta.columns else np.array(["NA"] * len(meta))
    fc = np.flatnonzero(groups == "FC")
    ref_mask = groups != "FC"

    n_components = min(30, z.shape[0] - 1, z.shape[1])
    pcs = PCA(n_components=n_components, random_state=0).fit_transform(z)

    # 1) k sweep on the full pool (one fit at max k, then prefix slices)
    ref = np.flatnonzero(ref_mask)
    nn = NearestNeighbors(n_neighbors=max(K_SWEEP)).fit(pcs[ref])
    nb = nn.kneighbors(pcs[fc], return_distance=False)
    base = 100 * float((groups[ref] == "BG").mean())
    for k in K_SWEEP:
        lab = groups[ref[nb[:, :k]]]
        rows.append(dict(dataset=key, variant=f"full_pool_k{k}", k=k,
                         frac_BG=100 * float((lab == "BG").mean()), base_BG=base,
                         n_FC_used=len(fc), n_FC_skipped=0))

    # 2) leave-own-section-out
    masks = [ref_mask & (sections != sections[i]) for i in fc]
    f, b, used, skip = frac_bg_restricted(pcs, groups, fc, masks, K_MAIN)
    rows.append(dict(dataset=key, variant="leave_own_section_out", k=K_MAIN,
                     frac_BG=f, base_BG=b, n_FC_used=used, n_FC_skipped=skip))

    # 3) same animal, other sections
    masks = [ref_mask & (animals == animals[i]) & (sections != sections[i]) for i in fc]
    f, b, used, skip = frac_bg_restricted(pcs, groups, fc, masks, K_MAIN)
    rows.append(dict(dataset=key, variant="same_animal_cross_section", k=K_MAIN,
                     frac_BG=f, base_BG=b, n_FC_used=used, n_FC_skipped=skip))

out = pd.DataFrame(rows)
out.to_csv(TABLES / "knn_sensitivity.csv", index=False)
print(out.to_string(index=False))
