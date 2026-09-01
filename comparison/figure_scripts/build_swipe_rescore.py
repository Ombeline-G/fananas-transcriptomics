#!/usr/bin/env python3
"""Second-rater re-scoring sample for the swipe validation (review point 26).

Draws a seed-fixed, stratified 15% sample of all swipe candidates (stratified
by dataset x original decision so both accepted and rejected cases are
represented), shuffles the order, and writes it WITHOUT the original
decisions. A second rater scores these candidates in the swipe tools; Cohen's
kappa is then computed with compute_rescore_kappa.py.

Output: tables/swipe_rescore_sample.csv (dataset, cell_label, section).
"""
import pandas as pd

from project_paths import ROOT, ZHUANG, HAO, SUYAMA

TABLES = ROOT / "comparison" / "tables"
SOURCES = {
    # MERFISH is two decision files since the ABCA-1 integration (2026-08-18): the
    # three originally validated series and the fourth animal. Both must be sampled,
    # otherwise the stated 15% covers only part of the candidate pool and the animal
    # carrying 172 of the 362 FC never enters the reliability check.
    "zhuang": ZHUANG / "swipe_validation" / "fc_decisions_zhuang.csv",
    "zhuang_abca1": ZHUANG / "swipe_validation" / "fc_decisions_abca1.csv",
    "hao": HAO / "swipe_validation" / "fc_decisions_hao.csv",
    "suyama": SUYAMA / "results" / "fc_decisions_suyama.csv",
}

parts = []
for key, path in SOURCES.items():
    d = pd.read_csv(path, dtype={"cell_label": str})
    take = (d.groupby("decision", group_keys=False)
             .apply(lambda g: g.sample(frac=0.15, random_state=0), include_groups=False))
    take = d.loc[take.index, ["cell_label", "section"]].copy()
    take.insert(0, "dataset", key)
    parts.append(take)
out = (pd.concat(parts, ignore_index=True)
         .sample(frac=1.0, random_state=1)      # shuffle presentation order
         .reset_index(drop=True))
out.to_csv(TABLES / "swipe_rescore_sample.csv", index=False)
print(out.groupby("dataset").size().to_dict(), "| total:", len(out))
