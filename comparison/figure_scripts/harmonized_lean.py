#!/usr/bin/env python3
"""Harmonized FC-vs-BG lean analysis — identical protocol across the three datasets.

Protocol pre-specified in comparison/notes/HARMONIZATION_protocol.md (2026-08-10):
purity gate (contam_Purkinje<1 & contam_granule<1, both groups) -> capture proxy
(n detected genes) -> 1:1 intra-section nearest-neighbour matching (10% relative
caliper, no replacement) -> paired Wilcoxon on per-dataset max-valid astro module
(dataset's own gate genes excluded) and on shared genes individually.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from project_paths import ROOT as PROJECT_ROOT
ROOT = str(PROJECT_ROOT)
TAB = os.path.join(ROOT, "comparison", "tables")

MODULES = {  # per-dataset max-valid astro module (gate genes of that dataset excluded)
    "zhuang": ["Aqp4", "Agt", "Gja1", "Aldh1l1"],
    "hao": ["Agt", "Slc1a2", "Slc1a3", "S100b", "Gjb6"],
    "suyama": ["Aqp4", "Agt", "Gja1", "Slc1a3", "S100b", "Aldh1l1", "Gjb6"],
}
SHARED_GENES = ["Aqp4", "S100b", "Gfap", "Ntsr2", "Aldh1l1", "Gja1", "Agt", "Stk32a"]
HAO_GATES = {"Aqp4", "Gja1", "Glul", "Aldh1l1"}
CALIPER = 0.10  # relative, on n detected genes


def match_pairs(fc, bg, capture):
    pairs = []
    for sec, f in fc.groupby("section"):
        b = bg[bg.section == sec]
        if b.empty:
            continue
        used = set()
        for i in f.index:
            d = (capture[b.index] - capture[i]).abs()
            d = d[[j for j in d.index if j not in used]]
            if len(d) == 0:
                continue
            j = d.idxmin()
            if d[j] <= CALIPER * max(capture[i], 1):
                used.add(j)
                pairs.append((i, j))
    return pairs


def main():
    mod_rows, gene_rows, persec_all, pairs_all = [], [], [], []
    for ds in ["zhuang", "hao", "suyama"]:
        c = pd.read_parquet(f"{TAB}/{ds}_astro_cells.parquet")
        e = pd.read_parquet(f"{TAB}/{ds}_astro_expr.parquet")
        if "index" in c.columns:
            c = c.set_index("index")
        e = e.set_index(c.index)

        z = (e - e.mean()) / e.std(ddof=0)  # z over ALL astroglia, pre-matching
        capture = (e > 0).sum(axis=1).astype(float)

        pure = c[(c.contam_Purkinje < 1) & (c.contam_granule < 1)]
        fc, bg = pure[pure.group == "FC"], pure[pure.group == "BG"]
        pairs = match_pairs(fc, bg, capture)
        fi = [p[0] for p in pairs]
        bi = [p[1] for p in pairs]

        mod = [g for g in MODULES[ds] if g in e.columns]
        score = z[mod].mean(axis=1)
        diffs = score[fi].values - score[bi].values
        p = wilcoxon(diffs).pvalue if len(diffs) >= 6 else np.nan
        sec_sign = pd.Series(diffs, index=fc.loc[fi, "section"].values).groupby(level=0).mean()
        # raw (unmatched, purity-gated only) delta for reference
        raw_delta = float(score[fc.index].mean() - score[bg.index].mean())
        mod_rows.append(dict(
            dataset=ds, n_FC_total=int((c.group == "FC").sum()), n_FC_pure=len(fc),
            n_pairs=len(pairs), module_genes=";".join(mod),
            capture_FC=float(capture[fi].median()), capture_BG=float(capture[bi].median()),
            delta_raw_pure=raw_delta,
            delta_module=float(np.mean(diffs)), p_wilcoxon=float(p),
            sections_pos=f"{int((sec_sign > 0).sum())}/{len(sec_sign)}",
        ))
        # per-section paired deltas (harmonized replacement for figS8_persection_lean.csv)
        persec = pd.DataFrame({
            "section": fc.loc[fi, "section"].values,
            "animal": fc.loc[fi, "animal"].values if "animal" in fc.columns else "NA",
            "pair_diff": diffs,
        }).groupby(["section", "animal"], as_index=False).agg(n_pairs=("pair_diff", "size"),
                                                             delta=("pair_diff", "mean"))
        persec.insert(0, "dataset", ds)
        persec_all.append(persec)
        # per-pair deltas (review point 16: cells -> sections -> animals hierarchy in Fig. 5a)
        pairs_all.append(pd.DataFrame({
            "dataset": ds,
            "section": fc.loc[fi, "section"].values,
            "animal": fc.loc[fi, "animal"].values if "animal" in fc.columns else "NA",
            "pair_diff": diffs,
        }))
        for g in SHARED_GENES:
            if g not in e.columns:
                continue
            gd = z[g][fi].values - z[g][bi].values
            gp = wilcoxon(gd).pvalue if len(gd) >= 6 and np.abs(gd).sum() > 0 else np.nan
            gene_rows.append(dict(
                dataset=ds, gene=g, delta_z=float(np.mean(gd)), p_wilcoxon=float(gp),
                interpretable=("NO (gate gene)" if (ds == "hao" and g in HAO_GATES) else "yes"),
            ))

    pd.DataFrame(mod_rows).to_csv(f"{TAB}/harmonized_lean_module.csv", index=False)
    pd.DataFrame(gene_rows).to_csv(f"{TAB}/harmonized_lean_genes.csv", index=False)
    pd.concat(persec_all).to_csv(f"{TAB}/harmonized_lean_persection.csv", index=False)
    pd.concat(pairs_all).to_csv(f"{TAB}/harmonized_lean_pairs.csv", index=False)
    print(pd.DataFrame(mod_rows).to_string(index=False))
    print()
    gp = pd.DataFrame(gene_rows).pivot(index="gene", columns="dataset", values="delta_z").round(3)
    pv = pd.DataFrame(gene_rows).pivot(index="gene", columns="dataset", values="p_wilcoxon")
    print("per-gene matched delta_z (FC-BG):")
    print(gp.to_string())
    print("\nper-gene p:")
    print(pv.round(4).to_string())


if __name__ == "__main__":
    main()
