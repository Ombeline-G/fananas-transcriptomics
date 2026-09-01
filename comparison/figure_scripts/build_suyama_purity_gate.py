#!/usr/bin/env python3
"""Build the canonical Suyama FC set (purity + anti-PCL) and capture-matched pairs.

The symmetric purity gate keeps cells with Purkinje and granule contamination
z-scores below 1 in both FC and BG. The anti-PCL rule then requires the FC soma
to sit at >= ANTI_PCL_UM (30 um) from the nearest Purkinje soma (d_pkj), so that
a slightly displaced BG cannot pass as FC. FC are matched 1:1 to BG by capture
(10% caliper). The script also exports the matched metabolic-pathway contrasts
used by Figure 5, plus a stricter 45-um control set.
"""
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from project_paths import ROOT as PROJECT_ROOT
ROOT = str(PROJECT_ROOT)
SU = os.path.join(ROOT, "datasets/suyama_2026_xenium", "results", "tables")
TAB = os.path.join(ROOT, "comparison", "tables")

CONTAM = {
    "contam_Purkinje": ["Pcp2", "Calb1", "Car8", "Itpr1", "Pcp4", "Prkcg"],
    "contam_granule": ["Gabra6", "Neurod1", "Rbfox3", "Cadps2", "Gng13"],
}
ASTRO_MOD = ["Aqp4", "Agt", "Gja1", "Slc1a3", "S100b", "Aldh1l1", "Gjb6"]
METAB = {
    "Glutamate-glutamine cycle": ["Slc1a3", "Slc1a2", "Glul", "Gls", "Slc38a1", "Slc38a3", "Got1", "Glud1"],
    "Glycolysis": ["Hk1", "Hk2", "Gpi1", "Pfkl", "Aldoa", "Aldoc", "Gapdh", "Pgk1", "Eno1", "Pkm", "Ldha"],
    "Lactate shuttle (MCT)": ["Ldha", "Ldhb", "Slc16a1", "Slc16a3", "Slc16a4"],
    "Glycogen metabolism": ["Gys1", "Pygb", "Pygl", "Ugp2", "Gbe1", "Agl", "Ppp1r3c"],
    "OXPHOS / mitochondria": ["Ndufa4", "Ndufs1", "Sdha", "Sdhb", "Uqcrc1", "Uqcrc2", "Cox4i1", "Cox5a", "Cox6c", "Atp5f1a", "Atp5f1b"],
    "Fatty-acid / lipid": ["Fabp7", "Apoe", "Cd36", "Acsbg1", "Elovl2", "Acsl3", "Cpt1a", "Fasn", "Scd2", "Dbi"],
    "Cholesterol synthesis": ["Hmgcr", "Hmgcs1", "Sqle", "Dhcr7", "Fdps", "Idi1", "Msmo1", "Insig1", "Lss"],
    "Glutathione / antioxidant": ["Gclc", "Gclm", "Gss", "Gpx1", "Gpx4", "Gsr", "Sod1", "Sod2", "Cat", "Prdx6"],
    "Iron / metal handling": ["Fth1", "Ftl1", "Mt1", "Mt2", "Cp", "Slc40a1", "Tfrc"],
    "TCA / ketone support": ["Bdh1", "Oxct1", "Acat1", "Idh1", "Idh2", "Aco2", "Cs", "Mdh1"],
}


def main():
    df = pd.read_csv(f"{SU}/bergmann_full_expr.csv")
    df["cell"] = df["cell"].astype(str)
    contamination_reference = set(
        str(x) for x in json.load(open(f"{SU}/suyama_fc_contaminated.json"))
    )
    meta = df[["cell", "group"]].copy()
    e = df.drop(columns=["cell", "group"])
    e.index = df["cell"]
    z = (e - e.mean()) / e.std(ddof=0)

    scores = pd.DataFrame(index=e.index)
    for name, genes in CONTAM.items():
        g = [x for x in genes if x in e.columns]
        scores[name] = z[g].mean(axis=1)
        print(f"{name}: genes used {g}")
    capture = (e > 0).sum(axis=1).astype(float)
    grp = pd.Series(meta.group.values, index=e.index).replace({"BG_PCL": "BG"})

    gate = (scores.contam_Purkinje < 1) & (scores.contam_granule < 1)
    fc_all = grp[grp == "FC"].index
    removed = set(fc_all[~gate[fc_all]])
    kept_purity = [c for c in fc_all if gate[c]]
    overlap = removed & contamination_reference
    print(f"\nFC total: {len(fc_all)} | removed by purity gate: {len(removed)}")
    print(f"removed cells supported by contamination reference: {len(overlap)}")

    # Regle anti-PCL canonique : soma a >= ANTI_PCL_UM du soma de Purkinje le plus
    # proche (d_pkj), pour qu'une BG au soma legerement remonte ne passe pas pour
    # une FC. 30 um ~ 1,5 diametre de soma de Purkinje au-dessus de la rangee.
    ANTI_PCL_UM = 30.0
    lam = pd.read_csv(f"{SU}/adult_bergmann_laminar_purkinje.csv")
    lam["cell"] = lam["cell"].astype(str)
    dpk = lam.set_index("cell")["d_pkj"]
    kept_fc = [c for c in kept_purity if c in dpk.index and dpk[c] >= ANTI_PCL_UM]
    print(f"purity-gated FC: {len(kept_purity)} | removed by anti-PCL rule (d_pkj < {ANTI_PCL_UM:.0f} um): "
          f"{len(kept_purity) - len(kept_fc)}")
    print(f"canonical FC retained: {len(kept_fc)}")

    bg_ok = [c for c in grp[grp == "BG"].index if gate[c]]
    used, pairs = set(), []
    for c in kept_fc:
        d = (capture[bg_ok] - capture[c]).abs()
        d = d[[j for j in d.index if j not in used]]
        if len(d) == 0:
            continue
        j = d.idxmin()
        if d[j] <= 0.10 * capture[c]:
            used.add(j)
            pairs.append((c, j))
    fi = [p[0] for p in pairs]
    bi = [p[1] for p in pairs]
    print(f"matched pairs: {len(pairs)} | capture FC {capture[fi].median():.0f} vs BG {capture[bi].median():.0f}")

    mod = [g for g in ASTRO_MOD if g in e.columns]
    sc = z[mod].mean(axis=1)
    diffs = sc[fi].values - sc[bi].values
    print(f"\nASTRO MODULE ({';'.join(mod)}):")
    print(f"  purity-gated matched delta = {np.mean(diffs):+.3f}  wilcoxon p = {wilcoxon(diffs).pvalue:.4g}")
    for g in ["Aqp4", "Gfap", "Ntsr2", "Agt", "Gja1", "S100b"]:
        gd = z[g][fi].values - z[g][bi].values
        print(f"  {g}: {np.mean(gd):+.3f}  p={wilcoxon(gd).pvalue:.3g}")

    rows = []
    for name, genes in METAB.items():
        g = [x for x in genes if x in e.columns]
        if len(g) < 3:
            continue
        ps = z[g].mean(axis=1)
        d = ps[fi].values - ps[bi].values
        rows.append(dict(pathway=name, n_genes=len(g), delta=float(np.mean(d)),
                         p_wilcoxon=float(wilcoxon(d).pvalue)))
    met = pd.DataFrame(rows)
    met.to_csv(f"{TAB}/fig5_metabolic_matched_suyama.csv", index=False)
    print("\nMETABOLIC (purity-gated matched pairs):")
    print(met.round(3).to_string(index=False))

    pd.DataFrame(dict(cell=fi, matched_bg=bi)).to_csv(
        f"{TAB}/suyama_purity_matched_pairs.csv", index=False
    )
    # Stable canonical FC set consumed by suyama_build.py.
    json.dump(sorted(kept_fc), open(f"{SU}/suyama_fc_gate.json", "w"))
    print(f"\nwrote canonical FC list -> suyama_fc_gate.json ({len(kept_fc)} cells)")

    # ---- Controle strict 45 um (robustesse ; seuil FC de Zhuang/Hao) ----
    MIN_DEPTH_UM = 45.0
    fc_strict = [c for c in kept_fc if dpk[c] >= MIN_DEPTH_UM]
    used_s, pairs_s = set(), []
    for c in fc_strict:
        d = (capture[bg_ok] - capture[c]).abs()
        d = d[[j for j in d.index if j not in used_s]]
        if len(d) == 0:
            continue
        j = d.idxmin()
        if d[j] <= 0.10 * capture[c]:
            used_s.add(j)
            pairs_s.append((c, j))
    fsi = [p[0] for p in pairs_s]
    bsi = [p[1] for p in pairs_s]
    print(f"\nSTRICT anti-PCL gate (d_pkj >= {MIN_DEPTH_UM:.0f} um): {len(fc_strict)} FC, {len(pairs_s)} matched pairs")
    rows_s = []
    diffs_s = sc[fsi].values - sc[bsi].values
    rows_s.append(dict(contrast="astro_module", n_pairs=len(pairs_s),
                       delta=float(np.mean(diffs_s)),
                       p_wilcoxon=float(wilcoxon(diffs_s).pvalue)))
    print(f"  astro module delta = {np.mean(diffs_s):+.3f}  p = {wilcoxon(diffs_s).pvalue:.4g}")
    for g in ["Aqp4", "Gfap", "Ntsr2", "Agt", "Gja1", "S100b"]:
        gd = z[g][fsi].values - z[g][bsi].values
        rows_s.append(dict(contrast=g, n_pairs=len(pairs_s),
                           delta=float(np.mean(gd)),
                           p_wilcoxon=float(wilcoxon(gd).pvalue)))
        print(f"  {g}: {np.mean(gd):+.3f}  p={wilcoxon(gd).pvalue:.3g}")
    pd.DataFrame(rows_s).to_csv(f"{TAB}/suyama_strict45_summary.csv", index=False)
    json.dump(sorted(fc_strict), open(f"{SU}/suyama_fc_gate_strict45.json", "w"))
    print(f"wrote strict FC list -> suyama_fc_gate_strict45.json ({len(fc_strict)} cells)"
          f" + suyama_strict45_summary.csv")


if __name__ == "__main__":
    main()
