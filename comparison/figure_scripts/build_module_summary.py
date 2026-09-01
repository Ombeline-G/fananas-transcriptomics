#!/usr/bin/env python3
"""Complete functional + metabolic module summary (review point 27).

Replaces the former Figures S6/S7 pair with one supplementary table and one
cross-platform heatmap (FigureS6_S8_S10.py). Two estimators, kept as produced
and labelled: functional modules = unmatched FC-vs-BG difference of module
z-scores with cell-level Mann-Whitney (fig7_functional_scores.csv, all three
datasets); metabolic pathways = capture-matched paired Wilcoxon
(fig5_metabolic_matched_{hao,suyama}.csv; the MERFISH panel lacks metabolic
genes). BH correction within each dataset x class. The cross-platform verdict
per module: 'concordant' when padj < 0.05 with the same sign in at least two
datasets that measure it; otherwise 'not concordant' (or 'single-platform' /
'not measurable elsewhere').

Output: tables/TableS_functional_metabolic_modules.csv
"""
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from project_paths import ROOT

TABLES = ROOT / "comparison" / "tables"
rows = []

func = pd.read_csv(TABLES / "fig7_functional_scores.csv")
for _, r in func.iterrows():
    rows.append(dict(module_class="functional", module=r.module, dataset=r.dataset,
                     n_genes=int(r.n_genes), delta=float(r.delta),
                     pval=float(r.pval), padj=float(r.padj),
                     estimator="unmatched FC-BG module z, cell-level Mann-Whitney"))

met_s = pd.read_csv(TABLES / "fig5_metabolic_matched_suyama.csv")
padj_s = multipletests(met_s.p_wilcoxon, method="fdr_bh")[1]
for (_, r), q in zip(met_s.iterrows(), padj_s):
    rows.append(dict(module_class="metabolic", module=r.pathway, dataset="Xenium",
                     n_genes=int(r.n_genes), delta=float(r.delta),
                     pval=float(r.p_wilcoxon), padj=float(q),
                     estimator="capture-matched pairs, paired Wilcoxon"))
met_h = pd.read_csv(TABLES / "fig5_metabolic_matched_hao.csv")
for _, r in met_h.iterrows():
    rows.append(dict(module_class="metabolic", module=r.pathway, dataset="Stereo-seq",
                     n_genes=int(r.n_genes), delta=float(r.mean_log2FC),
                     pval=np.nan, padj=np.nan,
                     estimator=f"capture-equalized pseudobulk mean log2FC "
                               f"(no pathway-level test; {int(r.n_sig_genes)} significant genes)"))
met_modules = sorted(set(met_s.pathway) | set(met_h.pathway))
for module in met_modules:
    rows.append(dict(module_class="metabolic", module=module, dataset="MERFISH",
                     n_genes=0, delta=np.nan, pval=np.nan, padj=np.nan,
                     estimator="not measurable (panel lacks metabolic genes)"))

out = pd.DataFrame(rows)

def verdict(sub):
    sig = sub[(sub.padj < 0.05) & sub.delta.notna()]
    measured = int(sub.delta.notna().sum())
    if len(sig) >= 2 and (np.sign(sig.delta).nunique() == 1):
        return f"concordant ({len(sig)}/{measured} datasets, same sign)"
    if len(sig) >= 2:
        return "significant but OPPOSITE signs"
    if len(sig) == 1:
        return f"single-platform ({sig.dataset.iloc[0].split(' ')[0]})"
    return "no significant dataset"

verdicts = out.groupby(["module_class", "module"], group_keys=False).apply(
    verdict, include_groups=False)
out["cross_platform_verdict"] = out.set_index(["module_class", "module"]).index.map(verdicts)
out.to_csv(TABLES / "TableS_functional_metabolic_modules.csv", index=False)

summary = out.drop_duplicates(["module_class", "module"])[
    ["module_class", "module", "cross_platform_verdict"]]
print(summary.to_string(index=False))
conc = summary[summary.cross_platform_verdict.str.startswith("concordant")]
print(f"\nconcordant modules: {len(conc)}")
if len(conc):
    for _, r in conc.iterrows():
        sub = out[(out.module == r.module) & (out.padj < 0.05)]
        print(r.module, "->", dict(zip(sub.dataset, sub.delta.round(3))))
