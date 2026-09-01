#!/usr/bin/env python3
"""Cohen's kappa between the original swipe decisions and a second rater's
re-scoring of tables/swipe_rescore_sample.csv (review point 26).

Usage: compute_rescore_kappa.py <rater2.csv>
where rater2.csv has columns: dataset, cell_label, decision (yes/no).
"""
import sys

import pandas as pd
from sklearn.metrics import cohen_kappa_score

from project_paths import ROOT, ZHUANG, HAO, SUYAMA

TABLES = ROOT / "comparison" / "tables"
SOURCES = {
    # Must mirror build_swipe_rescore.py: the sample is drawn from FOUR decision
    # files since the ABCA-1 integration (2026-08-18). Omitting zhuang_abca1 here
    # silently dropped its 61 sampled candidates from the merge (fixed 2026-08-20,
    # review Ombeline 3 point 13).
    "zhuang": ZHUANG / "swipe_validation" / "fc_decisions_zhuang.csv",
    "zhuang_abca1": ZHUANG / "swipe_validation" / "fc_decisions_abca1.csv",
    "hao": HAO / "swipe_validation" / "fc_decisions_hao.csv",
    "suyama": SUYAMA / "results" / "fc_decisions_suyama.csv",
}

rater2 = pd.read_csv(sys.argv[1], dtype={"cell_label": str})
sample = pd.read_csv(TABLES / "swipe_rescore_sample.csv", dtype={"cell_label": str})
orig = pd.concat([
    pd.read_csv(path, dtype={"cell_label": str})[["cell_label", "decision"]]
      .assign(dataset=key) for key, path in SOURCES.items()])

m1 = sample.merge(orig, on=["dataset", "cell_label"])
assert len(m1) == len(sample), (
    f"{len(sample) - len(m1)} sampled candidates have no original decision — "
    "SOURCES out of sync with build_swipe_rescore.py")
bad = set(rater2.decision.dropna()) - {"yes", "no"}
assert not bad, f"rater-2 decisions must resolve to yes/no before export, found: {bad}"
m = m1.merge(rater2, on=["dataset", "cell_label"], suffixes=("_r1", "_r2"))
missing = len(sample) - len(m)
if missing:
    print(f"WARNING: {missing} sampled candidates not scored by rater 2")
for key, sub in m.groupby("dataset"):
    k = cohen_kappa_score(sub.decision_r1, sub.decision_r2)
    agree = float((sub.decision_r1 == sub.decision_r2).mean())
    print(f"{key}: n={len(sub)}  agreement={agree:.2%}  Cohen kappa={k:.3f}")
k = cohen_kappa_score(m.decision_r1, m.decision_r2)
print(f"pooled: n={len(m)}  agreement={(m.decision_r1 == m.decision_r2).mean():.2%}  kappa={k:.3f}")
