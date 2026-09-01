#!/usr/bin/env python3
"""Restrict the active Hao position-first census table and summary to Mouse1."""
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from panels import module_scores
from project_paths import ROOT

TABLES = ROOT / "comparison" / "tables"
EXPR_PATH = TABLES / "hao_ML_census_expr.parquet"
SUMMARY_PATH = TABLES / "hao_ML_census_summary.csv"

expr = pd.read_parquet(EXPR_PATH)
mouse_column = "mouse" if "mouse" in expr.columns else "animal"
expr = expr[expr[mouse_column].astype(str).eq("Mouse1")].reset_index(drop=True)
assert len(expr) and set(expr[mouse_column].astype(str)) == {"Mouse1"}
expr.to_parquet(EXPR_PATH)

numeric = expr.select_dtypes(include="number").drop(
    columns=[column for column in ("nCount",) if column in expr.columns]
)
scores, _ = module_scores(numeric)
scores["group"] = expr.group.values
rows = []
for group, cells in expr.groupby("group"):
    row = {
        "group": group,
        "n": len(cells),
        "nCount_med": float(cells.nCount.median()),
        "berg_score": float(scores.loc[scores.group.eq(group), "bergmann"].median()),
        "astro_score": float(scores.loc[scores.group.eq(group), "astrocyte"].median()),
    }
    for gene in ["Kcnb2", "Kcnip3", "Kcnip1"]:
        row[f"{gene}_pct"] = float((cells[gene] > 0).mean() * 100)
        row[f"{gene}_mean"] = float(cells[gene].mean())
    rows.append(row)
summary = pd.DataFrame(rows)
summary.to_csv(SUMMARY_PATH, index=False)
print(summary.to_string(index=False))
