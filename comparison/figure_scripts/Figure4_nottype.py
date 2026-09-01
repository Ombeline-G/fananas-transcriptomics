#!/usr/bin/env python
"""Figure 4 — FC-BG divergence is bounded far below astrocyte-class divergence.

Panel order follows the narrative of Results §4 (review points 7-9):
a: per-cell pairwise silhouette distributions in the 30-PC kNN space\n   (violins + mean ticks; review points 8, 25).
b: MERFISH volcano, FC vs capture-matched BG (matched on n detected genes,
   same section — no spatial coordinate; review point 17) — the full landscape with the 9
   hits displayed and colored by predefined audit class, not removed
   (review point 9; fig4_DE_zhuang_depthmatched.csv + DE_audit_candidates.csv).
c: Xenium volcano, 44 matched pairs — no gene reaches padj < 0.05.
d: DE hit board — all 27 statistical hits NAMED as class-colored tiles
   (Zhuang 9 canonical; Hao 14 hatched, audit-excluded; Suyama 0; bold
   border = the Car8/Calb1 two-platform intersection; review points 9, 25).
e: equivalence magnitudes (review point 7) — ECDF of |Delta-log2| for the three
   contrasts on the MERFISH expressed-gene universe (fig4_equivalence_zhuang.csv);
   magnitudes are interpretable only on the single-soma multi-section platform
   (fig4_equivalence_platform_scoping.csv, Methods).
f: TOST caterpillar — the 324 per-gene 90% intervals ranked by effect,\n   colored by equivalence status (review points 7, 25).
The dissociated-atlas consistency check was demoted to Figure S12 (review
point 12): its FC-shift score can only detect the observed astrocytic shift,
so it is supplementary support, not a pillar.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, SEMANTIC, DSL, DS_ORDER, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"

from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize=(W2, 6.0))
gs = GridSpec(2, 12, figure=fig, wspace=3.2, hspace=0.55, top=0.935, bottom=0.085)
axA = fig.add_subplot(gs[0, 0:4])
axVz = fig.add_subplot(gs[0, 4:8]); axVs = fig.add_subplot(gs[0, 8:12])
axD = fig.add_subplot(gs[1, 0:5]); axEq = fig.add_subplot(gs[1, 5:9])
axT = fig.add_subplot(gs[1, 9:12])
# f carries a y-axis label and a 4-line annotation; nudge it right so it does not crowd e
_pt = axT.get_position()
axT.set_position([_pt.x0 + 0.030, _pt.y0, _pt.width - 0.030, _pt.height])

CLS_COLORS = {"astro_lean": (GROUP["FC"], "astro lean"),
              "contam_Purkinje": (SEMANTIC["contam"], "Purkinje contam."),
              "bergmann_program": (GROUP["BG"], "Bergmann"),
              "technical": ("#52514e", "technical"),
              "other": ("#b8b6b0", "other")}

# ---------- a: pairwise two-group silhouettes (review point 8) ----------
dist = pd.read_csv(TAB + "/fig4_distance.csv")
name_map = {"zhuang": "MERFISH", "hao": "Stereo-seq", "suyama": "Xenium"}
sil = {}
for pair_key in ["FC_BG", "BG_OtherCBastro", "FC_OtherCBastro"]:
    sub = dist[dist.pair == f"silhouette2_{pair_key}_pc30"].set_index("dataset").dist
    sil[pair_key] = [float(sub[name_map[k]]) for k in DS_ORDER]
assert np.allclose(sil["FC_BG"], [-0.0188354, 0.0698089, 0.0873872], atol=1e-4)
assert np.allclose(sil["BG_OtherCBastro"], [0.3088660, 0.1155596, 0.1467540], atol=1e-4)
assert np.allclose(sil["FC_OtherCBastro"], [0.1634790, 0.0706761, 0.1871799], atol=1e-4)
cells = pd.read_csv(TAB + "/silhouette_cells_pc30.csv")
x = np.arange(3); ws = 0.28
PAIRS = [("FC|BG", "#0b0b0b", "FC vs BG"),
         ("BG|OtherCBastro", GROUP["OtherCBastro"], "BG vs parenchymal"),
         ("FC|OtherCBastro", GROUP["FC"], "FC vs parenchymal")]
KEY2COL = {"FC|BG": "FC_BG", "BG|OtherCBastro": "BG_OtherCBastro",
           "FC|OtherCBastro": "FC_OtherCBastro"}
for off, (pair_key, color, lab) in zip([-ws, 0, ws], PAIRS):
    vals = [cells[(cells.dataset == ds) & (cells.pair == pair_key)].s_pc30.to_numpy()
            for ds in DS_ORDER]
    assert all(abs(v.mean() - sil[KEY2COL[pair_key]][i]) < 1e-4 for i, v in enumerate(vals))
    vp = axA.violinplot(vals, positions=x + off, widths=ws * 0.9, showextrema=False)
    for body in vp["bodies"]:
        body.set_facecolor(color); body.set_alpha(0.5); body.set_linewidth(0)
    for xi, v in zip(x + off, vals):
        axA.hlines(v.mean(), xi - 0.12, xi + 0.12, color=color, lw=2.0, zorder=5)
    axA.plot([], [], color=color, lw=3, alpha=0.6, label=lab)
axA.axhline(0, color="k", lw=0.8)
axA.axhspan(0.25, 0.92, color="#e1e0d9", alpha=0.5, lw=0)
axA.set_xticks(x); axA.set_xticklabels(["MERFISH", "Stereo-seq", "Xenium"], fontsize=7.5, rotation=90, va="top")
axA.set_ylabel("per-cell silhouette, two groups\n(30-PC space; tick = mean)")
axA.set_ylim(-0.78, 0.92)
axA.legend(frameon=False, fontsize=7, loc="lower right", labelspacing=0.28,
           handlelength=1.1, handletextpad=0.35, borderaxespad=0.05)
axA.set_title("Cluster separation", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axA, "a", dx=-0.18, dy=1.12)

# ---------- b, c: volcanoes, FC vs capture-matched BG (review points 9, 17) ----------
audit = pd.read_csv(TAB + "/DE_audit_candidates.csv")
zhits = audit[audit.dataset == "zhuang"].set_index("gene")
hhits = audit[audit.dataset == "hao"].set_index("gene")
assert len(zhits) == 9 and zhits.canonical_included.all()   # 4 animals: 13 -> 9 hits
assert len(hhits) == 14 and not hhits.canonical_included.any()
shared = sorted(set(zhits.index) & set(hhits.index))
assert shared == ["Calb1", "Car8"], shared

zv = pd.read_csv(TAB + "/fig4_DE_zhuang_depthmatched.csv")
sv = pd.read_csv(TAB + "/fig4_DE_suyama_depthmatched.csv")
assert int(((zv.padj < 0.05) & (zv.log2FC.abs() > 0.25)).sum()) == 9
assert int(((sv.padj < 0.05) & (sv.log2FC.abs() > 0.25)).sum()) == 0
assert float(sv.padj.min()) == 1.0
assert int((sv.pval < 0.05).sum()) == 145

# b: MERFISH — landscape + hits colored by audit class
yz = -np.log10(zv.pval.clip(lower=1e-15))
base = ~zv.gene.isin(zhits.index)
axVz.scatter(zv.log2FC[base], yz[base], s=4, c="#dddbd4", linewidths=0, zorder=2)
for cls, (col, lab) in CLS_COLORS.items():
    genes = zhits.index[zhits["class"] == cls]
    mask = zv.gene.isin(genes)
    if mask.any():
        axVz.scatter(zv.log2FC[mask], yz[mask], s=14, c=col, edgecolors="#333333",
                     linewidths=0.4, zorder=4, label=lab)
thr = float(zv[zv.padj < 0.05].pval.max())
axVz.axhline(-np.log10(thr), color="#898781", ls="--", lw=0.8)
axVz.text(-0.72, -np.log10(thr) + 0.25, "padj = 0.05", fontsize=7, color="#52514e")
for vx in (-0.25, 0.25):
    axVz.axvline(vx, color="#c3c2b7", lw=0.8, ls=":")
LABELS = {"Calb1": (0.05, -0.30, "left"), "Car8": (-0.04, -0.85, "right"),
          "Tnc": (-0.04, 0.55, "right"), "Gfap": (0.04, 0.35, "left"),
          "Fam107a": (-0.04, -0.45, "right"), "Ntsr2": (0.04, 0.40, "left"),
          "Agt": (-0.04, 0.10, "right"), "Aqp4": (0.04, -0.55, "left"),
          "Zeb2": (0.04, 0.10, "left")}   # the 9 hits (4 animals)
zv_idx = zv.set_index("gene")
for gname, (dx, dy, ha) in LABELS.items():
    gx = float(zv_idx.loc[gname, "log2FC"]); gy = float(-np.log10(zv_idx.loc[gname, "pval"]))
    axVz.text(gx + dx, gy + dy, f"$\\it{{{gname}}}$", fontsize=7, ha=ha, va="center")
axVz.set_xlim(-0.78, 0.78); axVz.set_ylim(-0.3, 14.2)
axVz.set_xlabel("Δlog2, FC vs matched BG")
axVz.set_ylabel("−log10 p (paired Wilcoxon)")
axVz.legend(frameon=False, fontsize=7, loc="upper right", labelspacing=0.25,
            handletextpad=0.25, borderaxespad=0.05, markerscale=1.1)
axVz.set_title("MERFISH — 335 pairs", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axVz, "b", dx=-0.22, dy=1.12)

# c: Xenium — flat landscape, nothing approaches significance
ys = -np.log10(sv.pval.clip(lower=1e-15))
axVs.scatter(sv.log2FC, ys, s=3, c="#dddbd4", linewidths=0, zorder=2)
for vx in (-0.25, 0.25):
    axVs.axvline(vx, color="#c3c2b7", lw=0.8, ls=":")
axVs.text(0.0, 12.4, "no gene reaches padj < 0.05\n(all padj = 1.0; 145/5,006\nnominal < 250 expected\nby chance)",
          ha="center", va="top", fontsize=7.2, color="#52514e")
axVs.set_xlim(-1.95, 1.95); axVs.set_ylim(-0.3, 14.2)
axVs.set_xlabel("Δlog2, FC vs matched BG")
axVs.set_ylabel("−log10 p (paired Wilcoxon)")
axVs.set_title("Xenium — 44 pairs", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axVs, "c", dx=-0.22, dy=1.12)

# ---------- d: DE hit board — every statistical hit named (review point 25) ----------
de = pd.read_csv(TAB + "/DE_canonical_summary.csv")
assert [int((de.dataset == key).sum()) for key in DS_ORDER] == [9, 0, 0]
counts = {"zhuang": zhits["class"].value_counts().to_dict(),
          "hao": hhits["class"].value_counts().to_dict()}
assert counts["zhuang"] == {"astro_lean": 4, "other": 2, "contam_Purkinje": 2, "bergmann_program": 1}
assert counts["hao"] == {"other": 6, "technical": 4, "contam_Purkinje": 4}
CLS_ORDER = ["astro_lean", "contam_Purkinje", "bergmann_program", "technical", "other"]
from matplotlib.patches import Rectangle
for col_i, (key, hits, hatched) in enumerate([("zhuang", zhits, False), ("hao", hhits, True)]):
    ordered = hits.reset_index().sort_values(
        ["class", "padj"], key=lambda c: c.map({k: i for i, k in enumerate(CLS_ORDER)})
        if c.name == "class" else c)
    for row_i, (_, r) in enumerate(ordered.iterrows()):
        y0 = 13.4 - row_i
        color = CLS_COLORS[r["class"]][0]
        shared_hit = r.gene in shared
        axD.add_patch(Rectangle((col_i * 1.12, y0), 1.0, 0.88, facecolor=color,
                                edgecolor="#0b0b0b" if shared_hit else "#333333",
                                lw=1.6 if shared_hit else 0.4,
                                hatch="///" if hatched else None, zorder=3))
        dark = r["class"] in ("astro_lean", "technical", "bergmann_program")
        axD.text(col_i * 1.12 + 0.5, y0 + 0.42, f"$\\it{{{r.gene}}}$", ha="center",
                 va="center", fontsize=7, zorder=4,
                 color="white" if dark else "#0b0b0b")
axD.text(0.5, 14.65, "MERFISH\n9 → 9", ha="center", va="bottom", fontsize=7.5)
axD.text(1.62, 14.65, "Stereo-seq\n14 → 0", ha="center", va="bottom", fontsize=7.5)
axD.text(2.85, 14.65, "Xenium\n0", ha="center", va="bottom", fontsize=7.5)
axD.text(2.28, 12.9, "colors: as in b\nhatched: audit-\nexcluded\nbold border:\nshared hits\n(all 3: none)",
         ha="left", va="top", fontsize=7, color="#0b0b0b")

axD.set_xlim(-0.15, 3.45); axD.set_ylim(-0.4, 17.4)
axD.set_xticks([]); axD.set_yticks([])
for sp in axD.spines.values():
    sp.set_visible(False)
axD.set_title("Every DE hit, named", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axD, "d", dx=-0.06, dy=1.06)

# ---------- e: effect-size distributions (equivalence magnitudes) ----------
eq = pd.read_csv(TAB + "/fig4_equivalence_zhuang.csv").set_index("gene")
assert len(eq) == 304
med = {c: float(eq[c].abs().median()) for c in ["d_FC_BG", "d_BG_par", "d_FC_par"]}
assert np.allclose([med["d_FC_BG"], med["d_BG_par"], med["d_FC_par"]],
                   [0.0308004, 0.1244108, 0.1279576], atol=1e-4)
CURVES = [("d_FC_BG", "#0b0b0b", "-", "FC vs BG"),
          ("d_BG_par", GROUP["OtherCBastro"], "-", "BG vs parenchymal"),
          ("d_FC_par", GROUP["FC"], "--", "FC vs parenchymal")]
for col, color, ls, lab in CURVES:
    v = np.sort(eq[col].abs().values)
    axEq.step(np.concatenate([[0], v]), np.linspace(0, 1, len(v) + 1),
              color=color, ls=ls, lw=1.6, where="post",
              label=f"{lab} ({med[col]:.2f})")
for vx in (0.25, 0.5):
    axEq.axvline(vx, color="#c3c2b7", lw=0.9, ls=":")
axEq.set_xlim(0, 1.2); axEq.set_ylim(0, 1.02)
axEq.set_xlabel("|Δlog2| per gene (n = 304)")
axEq.set_ylabel("cumulative fraction of genes")
axEq.legend(frameon=False, fontsize=7, loc="lower right", handlelength=1.2,
            labelspacing=0.4, bbox_to_anchor=(0.90, 0.02), title="median |Δlog2|",
            title_fontsize=7)
axEq.set_title("Magnitude of divergence", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axEq, "e", dx=-0.20, dy=1.12)

# ---------- f: TOST caterpillar — all 324 per-gene CIs (review point 25) ----------
tost = {}
for label in ["FC_BG", "BG_par"]:
    tost[label] = [int(eq[f"equiv_{m}_{label}"].sum()) for m in ("0.25", "0.5")]
assert tost["FC_BG"] == [277, 300] and tost["BG_par"] == [203, 250]
srt = eq.sort_values("sec_mean_FC_BG").reset_index()
xs = np.arange(len(srt))
in25 = srt["equiv_0.25_FC_BG"].to_numpy()
in50 = srt["equiv_0.5_FC_BG"].to_numpy()
col = np.where(in25, "#b8b6b0", np.where(in50, SEMANTIC["contam"], GROUP["FC"]))
axT.vlines(xs, srt.sec_lo90_FC_BG, srt.sec_hi90_FC_BG, colors=col, lw=0.7, zorder=3)
axT.scatter(xs, srt.sec_mean_FC_BG, s=1.5, c="#52514e", zorder=4)
for margin, ls in [(0.25, ":"), (0.5, "--")]:
    axT.axhline(margin, color="#898781", ls=ls, lw=0.8)
    axT.axhline(-margin, color="#898781", ls=ls, lw=0.8)
NUDGE_Y = {"Gfap": 0.10, "Aqp4": -0.08, "Agt": 0.0, "Calb1": 0.0, "Tnc": 0.02, "Car8": 0.04}
for _, r in srt[~in50].iterrows():
    xi = int(srt.index[srt.gene == r.gene][0])
    neg = r.sec_mean_FC_BG < 0
    axT.text(xi + (7 if neg else -7), r.sec_mean_FC_BG + NUDGE_Y.get(r.gene, 0),
             f"$\\it{{{r.gene}}}$", fontsize=7, ha="left" if neg else "right",
             va="center", color=GROUP["FC"])
axT.text(0.04, 0.97,
         "29 sections\n±0.25: 91% (277/304)\n±0.5: 99% (300/304)\nBG–parench.: 67% / 82%",
         transform=axT.transAxes, va="top", fontsize=7.2, color="#0b0b0b")
axT.set_xlabel("genes, ranked by effect", labelpad=6)
axT.set_ylabel("Δlog2, FC vs BG (90% CI)")
axT.set_xlim(-6, len(srt) + 5); axT.set_ylim(-1.05, 1.05)
axT.set_xticks([])
axT.set_title("Formal equivalence", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axT, "f", dx=-0.20, dy=1.12)

save(fig, "Figure4_nottype")
