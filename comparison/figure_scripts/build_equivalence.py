#!/usr/bin/env python3
"""Equivalence analysis of the FC-BG contrast (review point 7).

Absence of differential expression is not evidence of equivalence. This script
quantifies the MAGNITUDE of the FC-BG divergence against the bona fide
astrocyte-class divergence (BG vs parenchymal astrocytes), on the one platform
where per-cell magnitudes are interpretable: MERFISH (single-soma segmentation,
multi-section replication). Stereo-seq and Xenium magnitudes are excluded with
quantitative justification (build_equivalence_scoping.py): deconvolved-bin
effect sizes track detection breadth even after exact capture equalization, and
Xenium layer contrasts carry position-correlated neuropil spillover.

Definitions:
- universe: genes detected in >=10% of cells of at least one of the three
  groups (FC, BG, parenchymal astrocytes);
- effect size: pooled Delta-log2 = difference of group means of log2-normalized
  expression (the estimator used throughout the study). Sampling noise in the
  FC mean (n = 190) inflates |Delta|FC-BG, so the contrast with BG-parenchymal
  is conservative;
- TOST equivalence: per-section Delta (sections with >=3 FC, >=10 BG, >=10
  parenchymal — the same criteria as the centroid analysis), two one-sided
  tests at alpha 0.05 via the 90% t-interval across sections, margins
  |Delta-log2| < 0.25 (the DE effect threshold of this study) and < 0.5.
  Section-level intervals make equivalence HARDER to declare than the
  anti-conservatively narrow cell-level intervals (Statistics).

Output: tables/fig4_equivalence_zhuang.csv (one row per expressed gene).
"""
import numpy as np
import pandas as pd
from scipy import stats

from project_paths import ROOT

TABLES = ROOT / "comparison" / "tables"
GROUPS = ["FC", "BG", "OtherCBastro"]
MIN_DETECTION = 0.10
MIN_FC_SECTION, MIN_REF_SECTION = 3, 10
MARGINS = (0.25, 0.5)

meta = pd.read_parquet(TABLES / "zhuang_astro_cells.parquet")
expr = pd.read_parquet(TABLES / "zhuang_astro_expr.parquet").reindex(meta.cell)
groups = meta.group.to_numpy()
assert meta.group.value_counts().to_dict() == {"BG": 6084, "OtherCBastro": 1910, "FC": 362}   # 4 animals since ABCA-1 integration (2026-08-18)

detection = pd.DataFrame({g: (expr[groups == g] > 0).mean() for g in GROUPS})
expressed = detection.max(axis=1) >= MIN_DETECTION
E = expr.loc[:, expressed]
print(f"expressed-gene universe: {int(expressed.sum())}/{expr.shape[1]}")

mean = {g: E[groups == g].mean() for g in GROUPS}
out = pd.DataFrame({
    "gene": E.columns,
    "det_FC": detection.loc[expressed, "FC"].values,
    "det_BG": detection.loc[expressed, "BG"].values,
    "det_parenchymal": detection.loc[expressed, "OtherCBastro"].values,
    "d_FC_BG": (mean["FC"] - mean["BG"]).values,
    "d_BG_par": (mean["BG"] - mean["OtherCBastro"]).values,
    "d_FC_par": (mean["FC"] - mean["OtherCBastro"]).values,
}).set_index("gene")

# ---- section-level TOST ----
counts = meta.groupby("section").group.value_counts().unstack(fill_value=0)
sections = counts.index[(counts["FC"] >= MIN_FC_SECTION)
                        & (counts["BG"] >= MIN_REF_SECTION)
                        & (counts["OtherCBastro"] >= MIN_REF_SECTION)]
print(f"TOST sections: {len(sections)}")
section_arr = meta.section.to_numpy()
tcrit = stats.t.ppf(0.95, len(sections) - 1)
for label, (a, b) in [("FC_BG", ("FC", "BG")), ("BG_par", ("BG", "OtherCBastro"))]:
    per_sec = pd.DataFrame(
        [E[(section_arr == s) & (groups == a)].mean()
         - E[(section_arr == s) & (groups == b)].mean() for s in sections],
        index=sections)
    mu = per_sec.mean()
    se = per_sec.std(ddof=1) / np.sqrt(len(sections))
    out[f"sec_mean_{label}"] = mu
    out[f"sec_lo90_{label}"] = mu - tcrit * se
    out[f"sec_hi90_{label}"] = mu + tcrit * se
    for margin in MARGINS:
        eq = (out[f"sec_lo90_{label}"] > -margin) & (out[f"sec_hi90_{label}"] < margin)
        out[f"equiv_{margin}_{label}"] = eq
        print(f"{label}: TOST-equivalent within ±{margin}: "
              f"{int(eq.sum())}/{len(out)} ({100 * eq.mean():.1f}%)")

out.reset_index().to_csv(TABLES / "fig4_equivalence_zhuang.csv", index=False)

# ---- summary invariants ----
med = {c: float(out[c].abs().median()) for c in ["d_FC_BG", "d_BG_par", "d_FC_par"]}
frac25 = {c: float((out[c].abs() > 0.25).mean()) for c in ["d_FC_BG", "d_BG_par", "d_FC_par"]}
r_par = float(np.corrcoef(out.d_FC_par, out.d_BG_par)[0, 1])
not_eq_05 = sorted(out.index[~out["equiv_0.5_FC_BG"]])
print(f"median |d|: FC-BG {med['d_FC_BG']:.3f} · BG-par {med['d_BG_par']:.3f} · "
      f"FC-par {med['d_FC_par']:.3f}  (ratio {med['d_FC_BG'] / med['d_BG_par']:.2f})")
print(f"frac |d|>0.25: FC-BG {frac25['d_FC_BG']:.3f} · BG-par {frac25['d_BG_par']:.3f} · "
      f"FC-par {frac25['d_FC_par']:.3f}")
print(f"corr(d_FC_par, d_BG_par) over all expressed genes: r = {r_par:.3f}")
print(f"FC-BG non-equivalent at ±0.5: {not_eq_05}")
assert set(not_eq_05) <= {"Tnc", "Calb1", "Car8", "Fam107a"}, not_eq_05   # 4 animals: lean genes now bounded within ±0.5; Fam107a escapes
