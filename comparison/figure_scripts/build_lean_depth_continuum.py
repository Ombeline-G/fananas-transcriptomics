#!/usr/bin/env python3
"""Continuum test: astrocytic score versus soma position across ALL
Bergmann-class cells (review point 19).

Axis (MERFISH): signed perpendicular distance to the local Purkinje line, in
local nearest-neighbour units (scale-free; from the layer pipeline,
07_pool_layers_all4.py): positive = ML side (side called from the
MLI-versus-granule neighbourhood), negative = GCL side, 0 = the line.
Score: the dataset's max-valid astro module (Aqp4, Agt, Gja1, Aldh1l1),
z-scored over all astroglia — the same quantity as Figure 5.

Design constraints stated honestly:
- the census assigns every Bergmann-class ML soma beyond ~45 um to FC, so the
  zone between the PCL band and the threshold is sparsely populated and
  canonical BG barely exist on the ML side: the full-range per-section
  correlation is diluted by the sub-PCL mass (where the score is flat) and is
  reported as such;
- the informative summaries are the binned profile (per animal), the flat
  GCL side (negative control), and the within-FC gradient;
- sensitivity: rank-partial correlation controlling capture (n detected
  genes) and MLI/granule/Purkinje contamination scores.
Xenium (single section, unsigned distance, BG tightly at the PCL) is a
directional check only.

Outputs: tables/lean_depth_continuum_profile.csv (bins x animal),
         tables/lean_depth_continuum_sections.csv (per-section correlations).
"""
import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.stats import spearmanr

from project_paths import ROOT, ZHUANG

TABLES = ROOT / "comparison" / "tables"
MODULES = {"zhuang": ["Aqp4", "Agt", "Gja1", "Aldh1l1"],
           "suyama": ["Aqp4", "Agt", "Gja1", "Slc1a3", "S100b", "Aldh1l1", "Gjb6"]}
BINS = [-np.inf, -1, 1, 3, 6, 10, 15, np.inf]
BIN_LABELS = ["GCL<-1", "PCL-1..1", "1-3", "3-6", "6-10", "10-15", ">15"]


def load(key):
    m = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    E = pd.read_parquet(TABLES / f"{key}_astro_expr.parquet").reindex(m.cell)
    z = (E - E.mean()) / E.std(ddof=0)
    b = m[m.group.isin(["FC", "BG"])].copy()
    b["score"] = z.loc[:, MODULES[key]].mean(axis=1).values[m.group.isin(["FC", "BG"]).values]
    b["capture"] = (E > 0).sum(axis=1).values[m.group.isin(["FC", "BG"]).values]
    return b


def rank_partial_rho(a, bvals, covs):
    """Spearman of rank-residuals of a and bvals on the ranked covariates."""
    def resid(v):
        R = np.column_stack([np.ones(len(v))] + [pd.Series(c).rank().to_numpy() for c in covs])
        y = pd.Series(v).rank().to_numpy(dtype=float)
        beta, *_ = lstsq(R, y, rcond=None)
        return y - R @ beta
    return spearmanr(resid(a), resid(bvals))[0]


# ---------------- MERFISH ----------------
b = load("zhuang").set_index("cell")
cb = pd.read_parquet(ZHUANG / "cb_cells_layers_all4.parquet")
b = b.join(cb[["mli_frac", "gran_frac"]])
assert b.mli_frac.notna().all()
b["sdepth"] = np.where(b.mli_frac > b.gran_frac, 1.0, -1.0) * b.perp_norm
b["bin"] = pd.cut(b.sdepth, BINS, labels=BIN_LABELS)

profile_rows = []
for animal_label, frame in [("all", b)] + list(b.groupby("animal")):
    for bin_label, sub in frame.groupby("bin", observed=False):
        profile_rows.append(dict(
            dataset="zhuang", animal=animal_label, bin=bin_label, n=len(sub),
            n_FC=int((sub.group == "FC").sum()),
            mean_score=float(sub.score.mean()) if len(sub) else float("nan")))
profile = pd.DataFrame(profile_rows)

section_rows = []
for sec, sub in b.groupby("section"):
    if len(sub) < 30 or sub.sdepth.std() == 0:
        continue
    covs = [sub.capture, sub.contam_MLI, sub.contam_granule, sub.contam_Purkinje]
    fc = sub[sub.group == "FC"]
    section_rows.append(dict(
        dataset="zhuang", section=sec, animal=sub.animal.iloc[0], n=len(sub),
        n_FC=len(fc),
        rho_full=float(spearmanr(sub.score, sub.sdepth)[0]),
        rho_full_adjusted=float(rank_partial_rho(sub.score, sub.sdepth, covs)),
        rho_within_FC=float(spearmanr(fc.score, fc.sdepth)[0]) if len(fc) >= 8 else float("nan"),
    ))
sections = pd.DataFrame(section_rows)

profile.to_csv(TABLES / "lean_depth_continuum_profile.csv", index=False)
sections.to_csv(TABLES / "lean_depth_continuum_sections.csv", index=False)

allp = profile[profile.animal == "all"].set_index("bin").loc[BIN_LABELS]
print(allp[["n", "n_FC", "mean_score"]].round(3).to_string())
print(f"\nfull-range per-section rho: {int((sections.rho_full > 0).sum())}/{len(sections)} "
      f"positive, median {sections.rho_full.median():+.3f} "
      f"(adjusted {sections.rho_full_adjusted.median():+.3f}) — diluted by the sub-PCL mass")
fc = b[b.group == "FC"]
gcl = b[(b.group == "BG") & (b.sdepth < -1)]
pcl = b[b.bin == "PCL-1..1"]
print(f"GCL side mean {gcl.score.mean():+.3f} vs PCL band {pcl.score.mean():+.3f} (flat control)")
print(f"within FC: rho(score, sdepth) = {spearmanr(fc.score, fc.sdepth)[0]:+.3f} (n={len(fc)})")
for an, sub in fc.groupby("animal"):
    if len(sub) >= 8:
        print(f"  {an}: rho {spearmanr(sub.score, sub.sdepth)[0]:+.3f} (n={len(sub)})")
gap = b[(b.sdepth > 1) & (b.sdepth <= 6)]
print(f"sparse zone 1-6 units: n={len(gap)} ({int((gap.group=='FC').sum())} FC)")

# ---------------- Xenium (directional check) ----------------
s = load("suyama")
rho_all = spearmanr(s.score, s.depth_um)
sfc = s[s.group == "FC"]
print(f"\nXenium all Bergmann-class (n={len(s)}): rho(score, depth_um) = {rho_all[0]:+.3f}")
print(f"Xenium within FC (n={len(sfc)}): rho = {spearmanr(sfc.score, sfc.depth_um)[0]:+.3f}")
xb = pd.cut(s.depth_um, [0, 15, 30, 45, 70, np.inf], labels=["<15", "15-30", "30-45", "45-70", ">70"])
xprof = s.groupby(xb, observed=False).agg(n=("score", "size"),
                                          n_FC=("group", lambda g: int((g == "FC").sum())),
                                          mean_score=("score", "mean"))
xprof.insert(0, "animal", "all"); xprof.insert(0, "dataset", "suyama")
xprof = xprof.reset_index(names="bin")
pd.concat([profile, xprof]).to_csv(TABLES / "lean_depth_continuum_profile.csv", index=False)
print(xprof[["bin", "n", "n_FC", "mean_score"]].round(3).to_string(index=False))
