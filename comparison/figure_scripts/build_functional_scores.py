#!/usr/bin/env python3
"""Functional module scores, FC versus BG, on the three platforms (source of Figure S6 / Table S2).

For each dataset, every functional gene of func_panels.FUNC that the platform measures is
z-scored across the FC + BG cells of the canonical cohort; a module score is the mean z of
its measured genes; delta = mean(FC) - mean(BG); two-sided cell-level Mann-Whitney;
Benjamini-Hochberg within dataset.

Inputs (canonical cohorts, FC and BG only):
  tables/zhuang_func_expr.parquet, tables/suyama_func_expr.parquet  (zh_su_func_export.py)
  tables/hao_func_expr.parquet                                      (hao_func_export.py)
The MERFISH cohort must be the four-animal cohort of the study (asserted against
zhuang_astro_cells.parquet): until November 2026 the MERFISH rows of this table had been
left on the three-animal cohort that preceded the integration of ABCA-1.

Output: tables/fig7_functional_scores.csv (read by build_module_summary.py).
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from func_panels import FUNC
from project_paths import TABLES

DATASETS = [("zhuang", "MERFISH"), ("hao", "Stereo-seq"), ("suyama", "Xenium")]
META = ("cell", "index", "group", "animal")
rows = []
for key, name in DATASETS:
    f = pd.read_parquet(TABLES / f"{key}_func_expr.parquet")
    assert set(f.group) == {"FC", "BG"}, (key, set(f.group))
    # the functional export must hold the FC and BG of the canonical cohort: same FC, same
    # group sizes, and the same BG wherever BG is not subsampled (Stereo-seq BG is a per-section
    # random draw of the SURE PCL bins, made independently by hao_build.py and
    # hao_func_export.py: same population and size, different bins)
    cells = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    idc = "cell" if "cell" in cells.columns else "index"
    fid = "cell" if "cell" in f.columns else "index"
    ids = lambda t, col, g: set(t.loc[t.group == g, col].astype(str))
    assert ids(f, fid, "FC") == ids(cells, idc, "FC"), f"{key}: FC of the functional export != canonical FC"
    assert (f.group == "BG").sum() == (cells.group == "BG").sum(), f"{key}: BG count != canonical"
    if key != "hao":
        assert ids(f, fid, "BG") == ids(cells, idc, "BG"), f"{key}: BG of the functional export != canonical BG"
    genes = [g for g in f.columns if g not in META]
    z = (f[genes] - f[genes].mean()) / (f[genes].std() + 1e-9)
    fc, bg = (f.group == "FC").values, (f.group == "BG").values
    ds_rows = []
    for module, panel in FUNC.items():
        gg = [g for g in panel if g in genes]
        if not gg:
            continue                      # module not measured on this platform
        s = z[gg].mean(axis=1)
        ds_rows.append(dict(dataset=name, module=module,
                            delta=float(s[fc].mean() - s[bg].mean()),
                            mean_ML=float(s[fc].mean()), mean_PCL=float(s[bg].mean()),
                            pval=float(mannwhitneyu(s[fc], s[bg], alternative="two-sided").pvalue),
                            n_genes=len(gg)))
    padj = multipletests([r["pval"] for r in ds_rows], method="fdr_bh")[1]
    for r, q in zip(ds_rows, padj):
        r["padj"] = float(q)
    rows += ds_rows
    print(f"{name:10s} FC {int(fc.sum()):4d} · BG {int(bg.sum()):5d} · {len(ds_rows)}/{len(FUNC)} modules measured · "
          f"{int((padj < 0.05).sum())} with padj < 0.05")

out = pd.DataFrame(rows)[["dataset", "module", "delta", "mean_ML", "mean_PCL", "pval", "n_genes", "padj"]]
out.to_csv(TABLES / "fig7_functional_scores.csv", index=False)
print(out.round(4).to_string(index=False))
