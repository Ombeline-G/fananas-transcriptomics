#!/usr/bin/env python3
"""Rebuild the current full-panel and depth-matched FC-vs-BG DE tables.

Zhuang uses the same symmetric purity gate and within-section capture matching as
the harmonized lean analysis. Hao imports the current depth-equalized,
within-section permutation result. Suyama uses the canonical gate-56 pairs.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

from project_paths import ROOT as PROJECT_ROOT

ROOT = str(PROJECT_ROOT)
TAB = os.path.join(ROOT, "comparison", "tables")
HAO = os.path.join(ROOT, "datasets/hao_2024_stereo_seq")
SU = os.path.join(ROOT, "datasets/suyama_2026_xenium")
CALIPER = 0.10


def finalize(rows):
    out = pd.DataFrame(rows)
    out["padj"] = multipletests(out.pval.fillna(1), method="fdr_bh")[1]
    return out


def unpaired_de(expr, labels):
    a, b = expr[labels == "FC"], expr[labels == "BG"]
    rows = []
    for gene in expr.columns:
        x, y = a[gene].values, b[gene].values
        try:
            p = 1.0 if (x.std() == 0 and y.std() == 0) else mannwhitneyu(x, y).pvalue
        except Exception:
            p = 1.0
        rows.append(dict(gene=gene, log2FC=float(x.mean() - y.mean()), pval=p,
                         frac_ML=float((x > 0).mean()), frac_PCL=float((y > 0).mean()),
                         n_ML=len(x), n_PCL=len(y)))
    return finalize(rows)


def paired_de(expr, pairs):
    fi = [a for a, _ in pairs]
    bi = [b for _, b in pairs]
    rows = []
    for gene in expr.columns:
        x, y = expr.loc[fi, gene].values, expr.loc[bi, gene].values
        diff = x - y
        try:
            p = 1.0 if np.abs(diff).sum() == 0 else wilcoxon(diff).pvalue
        except Exception:
            p = 1.0
        rows.append(dict(gene=gene, log2FC=float(diff.mean()), pval=p,
                         frac_ML=float((x > 0).mean()), frac_PCL=float((y > 0).mean()),
                         n_ML=len(x), n_PCL=len(y)))
    return finalize(rows)


def capture_pairs(meta, expr):
    capture = (expr > 0).sum(axis=1).astype(float)
    pure = meta[(meta.contam_Purkinje < 1) & (meta.contam_granule < 1)]
    fc, bg = pure[pure.group == "FC"], pure[pure.group == "BG"]
    pairs = []
    for section, fsec in fc.groupby("section"):
        bsec = bg[bg.section == section]
        used = set()
        for cell in fsec.index:
            delta = (capture[bsec.index] - capture[cell]).abs()
            delta = delta[[i for i in delta.index if i not in used]]
            if delta.empty:
                continue
            partner = delta.idxmin()
            if delta[partner] <= CALIPER * max(capture[cell], 1):
                used.add(partner)
                pairs.append((cell, partner))
    return pairs


def n_sig(df):
    return int(((df.padj < 0.05) & (df.log2FC.abs() > 0.25)).sum())


# Zhuang: unmatched sensitivity + current harmonized 182-pair comparison.
zc = pd.read_parquet(os.path.join(TAB, "zhuang_astro_cells.parquet")).set_index("cell")
ze = pd.read_parquet(os.path.join(TAB, "zhuang_astro_expr.parquet")).reindex(zc.index)
zg = zc.group.values
zfull = unpaired_de(ze[zc.group.isin(["FC", "BG"])], zg[zc.group.isin(["FC", "BG"])])
zfull.to_csv(os.path.join(TAB, "fig4_DE_zhuang_full.csv"), index=False)
zpairs = capture_pairs(zc, ze)
zmatched = paired_de(ze, zpairs)
zmatched.to_csv(os.path.join(TAB, "fig4_DE_zhuang_depthmatched.csv"), index=False)
pd.DataFrame(zpairs, columns=["cell", "matched_bg"]).to_csv(
    os.path.join(TAB, "DE_zhuang_matched_pairs.csv"), index=False)
print(f"Zhuang: full {n_sig(zfull)} hits; matched {len(zpairs)} pairs, {n_sig(zmatched)} hits")

# Hao: current depth-matched result.
h = pd.read_csv(os.path.join(HAO, "analysis", "depth_matched", "DE_FC_vs_BG.csv"))
hm = pd.DataFrame({
    "gene": h.gene,
    "log2FC": h.log2FC,
    "pval": h.p_perm_section,
    "frac_ML": h.frac_ML,
    "frac_PCL": h.frac_PCL,
    "padj": h.padj,
    "n_ML": 76,
    "n_PCL": 1871,
})
hm.to_csv(os.path.join(TAB, "fig4_DE_hao_depthmatched.csv"), index=False)
print(f"Hao depth-matched: {n_sig(hm)} statistical hits; interpretation is handled by build_de_canonical.py")

# Suyama: unmatched gate-56 sensitivity + the 56 pre-specified capture pairs.
raw = pd.read_csv(os.path.join(SU, "results", "tables", "bergmann_full_expr.csv"))
raw["cell"] = raw.cell.astype(str)
raw["group"] = raw.group.replace({"BG_PCL": "BG"})
gate56 = set(map(str, json.load(open(os.path.join(SU, "results", "tables", "suyama_fc_gate.json")))))
raw = raw[(raw.group == "BG") | ((raw.group == "FC") & raw.cell.isin(gate56))].copy()
se = raw.set_index("cell").drop(columns="group")
se = se / np.log(2)  # Seurat ln(1+CP10k) -> log2(1+CP10k), exact — base homogenized 2026-08-17
sg = raw.set_index("cell").group.reindex(se.index).values
sfull = unpaired_de(se, sg)
sfull.to_csv(os.path.join(TAB, "fig4_DE_suyama_full.csv"), index=False)
sp = pd.read_csv(os.path.join(TAB, "suyama_purity_matched_pairs.csv"), dtype=str)
spairs = list(sp[["cell", "matched_bg"]].itertuples(index=False, name=None))
smatched = paired_de(se, spairs)
smatched.to_csv(os.path.join(TAB, "fig4_DE_suyama_depthmatched.csv"), index=False)
print(f"Suyama: full {n_sig(sfull)} hits; matched {len(spairs)} pairs, {n_sig(smatched)} hits")
