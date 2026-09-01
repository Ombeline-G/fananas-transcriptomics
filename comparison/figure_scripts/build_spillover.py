#!/usr/bin/env python3
"""Rebuild the spatial spillover-control table from the canonical cell parquets.

For each dataset, the nearest parenchymal-astrocyte distance is computed within
section. Coordinates are converted to micrometres before comparison. The final
OLS model is descriptive and adjusts the astrocyte score for distance plus the
five contamination scores already present in the canonical tables.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr, t as student_t

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_paths import ROOT


TABLES = ROOT / "comparison" / "tables"
DATASETS = {
    "zhuang": ("MERFISH", 1000.0),  # stored in millimetres
    "hao": ("Stereo-seq", 10.0),          # Stereo-seq bin coordinates
    "suyama": ("Xenium", 1.0),     # stored in micrometres
}
CONTAMINATION = [
    "contam_Purkinje",
    "contam_MLI",
    "contam_granule",
    "contam_oligoWM",
    "contam_vascular",
]


def ols_term(y: np.ndarray, predictors: pd.DataFrame, term: str) -> tuple[float, float]:
    """Return an OLS coefficient and two-sided t-test p-value without statsmodels."""
    columns = ["intercept", *predictors.columns]
    x = np.column_stack([np.ones(len(predictors)), predictors.to_numpy(dtype=float)])
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    dof = max(len(y) - x.shape[1], 1)
    sigma2 = float(residual @ residual) / dof
    covariance = sigma2 * np.linalg.pinv(x.T @ x)
    index = columns.index(term)
    se = float(np.sqrt(max(covariance[index, index], 0)))
    t_value = beta[index] / se if se > 0 else np.inf
    return float(beta[index]), float(2 * student_t.sf(abs(t_value), dof))


def add_nearest_parenchymal_distance(cells: pd.DataFrame, scale: float) -> pd.DataFrame:
    cells = cells.copy()
    cells["dpar_um"] = np.nan
    for _, section in cells.groupby("section"):
        reference = section.loc[section.group == "OtherCBastro", ["x", "y"]].to_numpy() * scale
        if not len(reference):
            continue
        query = section[["x", "y"]].to_numpy() * scale
        cells.loc[section.index, "dpar_um"] = cKDTree(reference).query(query, k=1)[0]
    return cells


def summarize(key: str, label: str, scale: float) -> dict:
    cells = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    cells = add_nearest_parenchymal_distance(cells, scale)
    fc = cells[(cells.group == "FC") & cells.dpar_um.notna()].copy()
    bg = cells[(cells.group == "BG") & cells.dpar_um.notna()].copy()

    rho, p_rho = spearmanr(fc["astrocyte"], fc["dpar_um"])
    model_data = pd.concat([fc, bg], ignore_index=True)
    model_data["is_FC"] = (model_data.group == "FC").astype(int)
    predictors = ["is_FC", "dpar_um", *CONTAMINATION]
    adjusted_delta, adjusted_p = ols_term(
        model_data["astrocyte"].to_numpy(dtype=float),
        model_data[predictors],
        "is_FC",
    )

    return {
        "dataset": label,
        "n_FC": len(fc),
        "dpar_FC_um": np.median(fc.dpar_um),
        "dpar_BG_um": np.median(bg.dpar_um),
        "p_parench_closer": mannwhitneyu(fc.dpar_um, bg.dpar_um).pvalue,
        "spearman_astro_vs_dist": rho,
        "p_prox": p_rho,
        "astro_delta_raw": fc.astrocyte.mean() - bg.astrocyte.mean(),
        "astro_delta_dpar_contam_adj": adjusted_delta,
        "p_adj": adjusted_p,
    }


rows = [summarize(key, label, scale) for key, (label, scale) in DATASETS.items()]
out = pd.DataFrame(rows)
out.to_csv(TABLES / "fig4S3_spillover.csv", index=False)

# per-section medians of the FC/BG distances to the nearest parenchymal astrocyte
# (review point 25: Figure 5c shows the underlying per-section data, not bars)
persec = []
for key, (label, scale) in DATASETS.items():
    cells = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    cells = add_nearest_parenchymal_distance(cells, scale)
    sub = cells[cells.group.isin(["FC", "BG"]) & cells.dpar_um.notna()]
    for sec, g in sub.groupby("section"):
        fc, bg = g[g.group == "FC"], g[g.group == "BG"]
        if len(fc) == 0 or len(bg) == 0:
            continue
        persec.append(dict(dataset=key, section=sec,
                           median_FC_um=float(fc.dpar_um.median()),
                           median_BG_um=float(bg.dpar_um.median()),
                           n_FC=len(fc), n_BG=len(bg)))
pd.DataFrame(persec).to_csv(TABLES / "spillover_persection.csv", index=False)
print(pd.DataFrame(persec).groupby("dataset").size().to_dict())
print(out.to_string(index=False))
