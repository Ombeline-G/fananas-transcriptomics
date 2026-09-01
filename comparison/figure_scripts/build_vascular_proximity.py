#!/usr/bin/env python3
"""Vascular proximity of FC versus BG, and Aqp4 versus vessel distance
(review point 21).

Vessel proxies: MERFISH — endothelial somata of the Allen taxonomy (subclass
'333 Endo NN'), per section, planar distance (mm converted to um); Xenium —
the four Vascular/stroma clusters of the adult section (c8/c20/c21/c23; a
mixed vascular-leptomeningeal proxy, stated as such). Both are 2D single-soma
proxies: they underestimate 3D vessel proximity and mark endothelial somata
rather than the vessel wall — stated in Methods.

Outputs: tables/vascular_proximity.csv
 - level='section' rows (MERFISH): median distance to the nearest endothelial
   soma for FC and BG per section;
 - level='summary' rows: pooled medians, the within-group Spearman
   correlations of Aqp4 (z) and of the astro-module score with vessel
   distance (cell-level, descriptive — Statistics), the cell-level
   Mann-Whitney p for FC-vs-BG distance, and (MERFISH) the section-level
   paired Wilcoxon p on per-section medians — the appropriate level given
   nesting. The two single-soma platforms disagree in DIRECTION, which is the
   reported result; neither supports a perivascular FC position.
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, mannwhitneyu, wilcoxon

from project_paths import ROOT, ZHUANG, SUYAMA

TABLES = ROOT / "comparison" / "tables"
MODULES = {"zhuang": ["Aqp4", "Agt", "Gja1", "Aldh1l1"],
           "suyama": ["Aqp4", "Agt", "Gja1", "Slc1a3", "S100b", "Aldh1l1", "Gjb6"]}
rows = []


def load(key):
    m = pd.read_parquet(TABLES / f"{key}_astro_cells.parquet")
    if "cell" not in m.columns:
        m = m.rename(columns={"index": "cell"})
    m["cell"] = m.cell.astype(str)
    E = pd.read_parquet(TABLES / f"{key}_astro_expr.parquet").reindex(m.cell)
    z = (E - E.mean()) / E.std(ddof=0)
    keep = m.group.isin(["FC", "BG"]).values
    b = m[keep].copy()
    b["aqp4"] = z["Aqp4"].values[keep]
    b["score"] = z[MODULES[key]].mean(axis=1).values[keep]
    return b


def summarize(key, b):
    fc, bg = b[b.group == "FC"], b[b.group == "BG"]
    rows.append(dict(
        dataset=key, level="summary", section="all",
        median_FC_um=float(fc.dvasc.median()), median_BG_um=float(bg.dvasc.median()),
        n_FC=len(fc), n_BG=len(bg),
        rho_aqp4_dvasc_FC=float(spearmanr(fc.aqp4, fc.dvasc)[0]),
        rho_score_dvasc_FC=float(spearmanr(fc.score, fc.dvasc)[0]),
        rho_aqp4_dvasc_BG=float(spearmanr(bg.aqp4, bg.dvasc)[0]),
        p_MW_cell=float(mannwhitneyu(fc.dvasc, bg.dvasc, alternative="two-sided")[1]),
    ))


# ---------------- MERFISH ----------------
b = load("zhuang")
cb = pd.read_parquet(ZHUANG / "cb_cells_layers_all4.parquet")
endo = cb[cb.subclass == "333 Endo NN"]
dists = np.full(len(b), np.nan)
for sec, sub in b.groupby("section"):
    ev = endo[endo.brain_section_label == sec]
    assert len(ev) >= 10, (sec, len(ev))
    tree = cKDTree(ev[["x", "y"]].to_numpy())
    d, _ = tree.query(sub[["x", "y"]].to_numpy())
    dists[b.index.get_indexer(sub.index)] = d * 1000.0  # mm -> um
b["dvasc"] = dists
assert b.dvasc.notna().all()
for sec, sub in b.groupby("section"):
    fc, bg = sub[sub.group == "FC"], sub[sub.group == "BG"]
    if len(fc) == 0:
        continue
    rows.append(dict(dataset="zhuang", level="section", section=sec,
                     median_FC_um=float(fc.dvasc.median()),
                     median_BG_um=float(bg.dvasc.median()),
                     n_FC=len(fc), n_BG=len(bg)))
summarize("zhuang", b)

# ---------------- Xenium ----------------
s = load("suyama")
sc = pd.read_csv(SUYAMA / "results" / "tables" / "adult_cells_scores.csv", dtype={"cell": str})
vasc = sc[sc.cl.isin(["c8", "c20", "c21", "c23"])]
tree = cKDTree(vasc[["x", "y"]].to_numpy())
d, _ = tree.query(s[["x", "y"]].to_numpy())
s["dvasc"] = d
summarize("suyama", s)

out = pd.DataFrame(rows)
# section-level paired test (the appropriate level for MERFISH; Statistics)
_sz = out[(out.level == "section") & (out.dataset == "zhuang")]
p_sec = float(wilcoxon(_sz.median_FC_um.values, _sz.median_BG_um.values)[1])
out.loc[out.dataset.eq("zhuang") & out.level.eq("summary"), "p_section_paired"] = p_sec
out.to_csv(TABLES / "vascular_proximity.csv", index=False)
summ = out[out.level == "summary"].set_index("dataset")
print(summ.round(3).to_string())
sec = out[(out.level == "section") & (out.dataset == "zhuang")]
closer = int((sec.median_FC_um < sec.median_BG_um).sum())
print(f"\nMERFISH per-section: FC median closer in {closer}/{len(sec)} sections"
      f" | paired Wilcoxon on section medians p = {p_sec:.3g}")
