#!/usr/bin/env python
"""Figure 2 — molecular-layer astroglial census (redesigned 2026-08-17, review point 3).

Panels: (a) full census map; (b) zoom on a WM boundary with Astro-NT chains;
(c) soma-to-WM distance; (d) composition raw -> ML proper (>30 µm from WM AND pia);
(e) transcriptomic identity of every census group (bergmann vs astrocyte module z);
(f) Bergmann share of ML astroglia across the three datasets.
The audit block asserts every displayed value and aborts on mismatch.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, ENTITY, panel_letter, save, W2
from geom_layers import layer_regions, paint_layers

use_style()
from project_paths import ROOT as PROJECT_ROOT
BASE = str(PROJECT_ROOT)
TAB = BASE + "/comparison/tables"
ZH = BASE + "/datasets/zhang_zhuang_2023_merfish"
SU = BASE + "/datasets/suyama_2026_xenium"

MAP_SEC = "Zhuang-ABCA-3.005"
C = {"Bergmann": GROUP["FC"], "Astro-NT": ENTITY["Astro_NT"], "Astro-CB": ENTITY["Astro_CB"]}
NT_EDGE = "#7a5200"
LAYER_TINT = {"ML": "#dcebf9", "PCL": "#f6d9a8", "GCL": "#e7ddf3", "WM": "#e4e7eb"}
PCL_LINE = "#b8560f"
ORDER3 = ["Bergmann", "Astro-NT", "Astro-CB"]
XLAB3 = ["displaced\nBergmann", "WM\nastro", "parench.\nastro"]

# ================= DATA AUDIT (aborts on mismatch) =================
print("== DATA AUDIT (Figure 2, redesigned) ==")
zc = pd.read_parquet(TAB + "/zhuang_ML_census_cells.parquet")
z45 = zc[zc.depth_from_PCL >= 45].copy()
z45["dist_border"] = z45[["dist_to_WM", "dist_to_pia"]].min(axis=1)
counts = {t: int((z45.type == t).sum()) for t in ORDER3}
assert counts == {"Bergmann": 263, "Astro-NT": 246, "Astro-CB": 10}, counts
TOT = sum(counts.values()); assert TOT == 519
meds = {t: float(z45.loc[z45.type == t, "dist_to_WM"].median()) for t in ORDER3}
assert round(meds["Bergmann"]) == 96 and round(meds["Astro-NT"]) == 18 and round(meds["Astro-CB"]) == 66
prop = z45[z45.dist_border > 30]
pcounts = {t: int((prop.type == t).sum()) for t in ORDER3}
assert pcounts == {"Bergmann": 214, "Astro-NT": 52, "Astro-CB": 6}, pcounts
PTOT = sum(pcounts.values()); assert PTOT == 272
ident = pd.read_csv(TAB + "/zhuang_census_identity_cells.csv")
isum = pd.read_csv(TAB + "/zhuang_census_identity.csv").set_index("group")
assert round(isum.loc["AstroNT_interior", "berg_z_med"], 2) == -0.46
assert round(isum.loc["AstroNT_interior", "astro_z_med"], 2) == 0.58
assert round(isum.loc["census_Bergmann", "berg_z_med"], 2) == 0.31
# cross-dataset shares
hc = pd.read_parquet(TAB + "/hao_ML_census_expr.parquet")
hgv = hc.group.value_counts().to_dict()
mlbg = int(hgv["FC"])              # ML Bergmann-class bins (census label "FC")
mlnb = int(hgv["Astro_ML_nonBG"])  # ML bins not assigned to the Bergmann class
assert (mlbg, mlnb) == (450, 114), hgv
sdec = pd.read_csv(SU + "/results/fc_decisions_suyama.csv")
n_swipe = int((sdec.decision == "yes").sum()); assert n_swipe == 73, n_swipe
# Explicit test of the geometric ML band (review point 3e / Ombeline point 1): the swipe pool
# was NOT a class-blind census — 150 Astrocyte-cluster somata sit in the band and never entered
# it. Panel f therefore shows the band -> folial two-stage share, with the folial denominator
# defined by that explicit test (0/150 folial: marginal chain n=116 + deep-nuclei interface n=34).
a150 = pd.read_csv(TAB + "/suyama_ML_astro150.csv")
assert len(a150) == 150 and a150.cl.value_counts().to_dict() == {"c11": 116, "c13": 34}
assert not set(a150.cell) & set(sdec.cell_label)          # 0/150 in the swipe candidate pool
SU_BAND = (n_swipe, n_swipe + len(a150))                  # 73 / 223 geometric-band somata
Z_RAW = 100 * counts["Bergmann"] / TOT                    # raw ML mask, Zhuang
Z_PRO = 100 * pcounts["Bergmann"] / PTOT                  # ML proper, Zhuang
H_BAND = 100 * mlbg / (mlbg + mlnb)
assert (round(Z_RAW, 1), round(Z_PRO, 1), round(H_BAND, 1)) == (50.7, 78.7, 79.8)
assert round(100 * SU_BAND[0] / SU_BAND[1], 1) == 32.7
print(f"  raw 263/246/10 ✓  ML-proper 214/52/6 ✓  identity medians ✓  Hao 450/114 ✓  "
      f"Suyama band 73/223 -> folial 73/73 (150 tested, 0 folial) ✓")

# map region
zz = pd.read_parquet(ZH + "/cb_cells_multivote.parquet")
zz = zz[zz.brain_section_label == MAP_SEC]
reg = layer_regions(zz.x, zz.y, zz.is_gran, zz.subclass.astype(str).str.contains("Oligo"), umpp=1000.0)
XX, YY, T = reg["XX"], reg["YY"], np.asarray(reg["tissue"])
gmap = z45[z45.section == MAP_SEC]

# ================= figure =================
fig = plt.figure(figsize=(W2, 6.1))
gs = GridSpec(2, 3, figure=fig, width_ratios=[1.12, 1.0, 1.05], height_ratios=[1.0, 0.92],
              wspace=0.50, hspace=0.42, top=0.945, bottom=0.165, left=0.075, right=0.985)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[0, 2])
axD = fig.add_subplot(gs[1, 0]); axE = fig.add_subplot(gs[1, 1]); axF = fig.add_subplot(gs[1, 2])


nt_all = zz[zz.subclass.astype(str).str.contains("Astro-NT")]
nt_ghost = nt_all[~nt_all.index.astype(str).isin(set(gmap.cell.astype(str)))]
print(f"  ghost Astro-NT outside the census on {MAP_SEC}: {len(nt_ghost)}")


def draw_census_map(ax, xlim=None, ylim=None, star=40, dot=16, ghost=6):
    paint_layers(ax, reg, colors=LAYER_TINT, pcl_line=PCL_LINE, pcl_lw=0.9)
    Ts = ndi.gaussian_filter(T.astype(float), sigma=1.6)
    ax.contour(XX, YY, Ts, levels=[0.5], colors=["#b6bdc4"], linewidths=0.55, zorder=8)
    # subclass home: every Astro-NT of the section OUTSIDE the ML-mask census (mostly WM/GCL)
    ax.scatter(nt_ghost.x, nt_ghost.y, s=ghost, facecolors="none", edgecolors="#c9962e",
               linewidths=0.5, alpha=0.6, zorder=4)
    for t in ["Astro-NT", "Astro-CB"]:
        gg = gmap[gmap.type == t]
        ec = NT_EDGE if t == "Astro-NT" else "#333333"
        ax.scatter(gg.x, gg.y, s=dot, c=C[t], edgecolors=ec, linewidths=0.45, zorder=6)
    gg = gmap[gmap.type == "Bergmann"]
    ax.scatter(gg.x, gg.y, s=star * 2, marker="*", c="white", edgecolors="none", zorder=5)
    ax.scatter(gg.x, gg.y, s=star, marker="*", c=C["Bergmann"], edgecolors="#7a1213",
               linewidths=0.45, alpha=0.95, zorder=6)
    if xlim: ax.set_xlim(*xlim)
    if ylim: ax.set_ylim(*ylim)
    ax.set_aspect("equal"); ax.axis("off")


# ---------- a: full census map ----------
pad = 0.15
x0, x1 = gmap.x.min() - pad, gmap.x.max() + pad
y0, y1 = gmap.y.min() - pad, gmap.y.max() + pad
cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
half = max(x1 - x0, y1 - y0) / 2
draw_census_map(axA, xlim=(cx - half, cx + half), ylim=(cy - half, cy + half))
xl0, xl1 = axA.get_xlim(); yl0, yl1 = sorted(axA.get_ylim())
L = 0.5
xs = xl1 - (xl1 - xl0) * 0.05 - L; ys = yl0 + (yl1 - yl0) * 0.05
axA.add_line(Line2D([xs, xs + L], [ys, ys], color="#0b0b0b", lw=1.8, solid_capstyle="butt", zorder=11))
from matplotlib import patheffects as pe
axA.text(xs + L / 2, ys + (yl1 - yl0) * 0.022, "500 µm", ha="center", va="bottom", fontsize=7,
         zorder=11, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
axA.set_title(MAP_SEC.replace("Zhuang-", ""), loc="right", fontsize=7.0, color="#5a5a57", pad=2)
panel_letter(axA, "a", dx=0.0, dy=1.04)

# zoom rectangle on a
BXL, BYL = (12.28, 12.92), (3.80, 4.34)
from matplotlib.patches import Rectangle
axA.add_patch(Rectangle((BXL[0], BYL[0]), BXL[1] - BXL[0], BYL[1] - BYL[0],
                        fill=False, edgecolor="#0b0b0b", lw=0.9, zorder=12))

# ---------- b: zoom on the WM boundary ----------
draw_census_map(axB, xlim=BXL, ylim=BYL, star=70, dot=34, ghost=14)
Lz = 0.1
xsz = BXL[1] - 0.03 - Lz; ysz = BYL[0] + 0.035
axB.add_line(Line2D([xsz, xsz + Lz], [ysz, ysz], color="#0b0b0b", lw=1.8, solid_capstyle="butt", zorder=11))
axB.text(xsz + Lz / 2, ysz + 0.012, "100 µm", ha="center", va="bottom", fontsize=7,
         zorder=11, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
axB.set_title("Astro-NT hug the deep tissue border;\nFC lie within the ML", fontsize=7.5, loc="right", pad=3)
panel_letter(axB, "b", dx=0.0, dy=1.04)

# ---------- c: soma-to-WM distance ----------
rng = np.random.default_rng(0)
dmax = 0
for i, t in enumerate(ORDER3):
    d = z45.loc[z45.type == t, "dist_to_WM"].values
    dmax = max(dmax, d.max())
    if len(d) >= 25:
        pc = axC.violinplot([d], positions=[i], showextrema=False, widths=0.75)["bodies"][0]
        pc.set_facecolor(C[t]); pc.set_alpha(0.15)
    fill = "#d08c00" if t == "Astro-NT" else C[t]
    edge = NT_EDGE if t == "Astro-NT" else "none"
    axC.scatter(i + rng.normal(0, 0.05, len(d)), d, s=5, c=fill, alpha=0.5, edgecolors=edge, linewidths=0.3)
    axC.hlines(np.median(d), i - 0.3, i + 0.3, color="#0b0b0b", lw=1.3, zorder=5)
axC.axhline(30, ls=(0, (4, 3)), c="#898781", lw=0.9)
axC.text(2.55, 34, "30 µm", fontsize=7.0, color="#898781", ha="right")
axC.set_xticks(range(3)); axC.set_xticklabels(XLAB3, fontsize=7.5)
axC.tick_params(axis="x", length=0)
axC.set_xlim(-0.6, 2.6); axC.set_ylim(0, min(dmax * 1.05, 430))
axC.set_ylabel("soma-to-white-matter\ndistance (µm)")
axC.set_yticks([0, 100, 200, 300, 400])
panel_letter(axC, "c", dx=-0.30, dy=1.06)

# ---------- d: composition raw -> ML proper ----------
w = 0.38
xs_ = np.arange(3)
raw_pct = [100 * counts[t] / TOT for t in ORDER3]
pro_pct = [100 * pcounts[t] / PTOT for t in ORDER3]
b1 = axD.bar(xs_ - w / 2, raw_pct, w, color=[C[t] for t in ORDER3], alpha=0.45,
             edgecolor="#555555", linewidth=0.5, label=f"raw ML mask (n={TOT})")
b2 = axD.bar(xs_ + w / 2, pro_pct, w, color=[C[t] for t in ORDER3], alpha=0.95,
             edgecolor="#222222", linewidth=0.7,
             label=f"ML proper: >30 µm from WM & pia (n={PTOT})")
for x, v, n in zip(xs_ - w / 2, raw_pct, [counts[t] for t in ORDER3]):
    axD.text(x, v + 1.2, f"{v:.0f}%", ha="center", fontsize=7.0, color="#555555")
for x, v, n in zip(xs_ + w / 2, pro_pct, [pcounts[t] for t in ORDER3]):
    axD.text(x, v + 1.2, f"{v:.0f}%", ha="center", fontsize=7.0, fontweight="bold")
axD.set_xticks(xs_); axD.set_xticklabels(XLAB3, fontsize=7.5)
axD.tick_params(axis="x", length=0)
axD.set_ylabel("share of ML-mask astroglia (%)")
axD.set_ylim(0, 100)
axD.legend(fontsize=7.0, frameon=False, loc="upper right", handlelength=1.2, borderaxespad=0.1)
panel_letter(axD, "d", dx=-0.30, dy=1.06)

# ---------- e: identity of the census groups ----------
zone_face = {"border": "#eda100", "interior": "#eda100"}
sub = ident[ident.type == "Astro-NT"]
for zone, ec, lw, lab in [("border", NT_EDGE, 0.4, None), ("interior", "#0b0b0b", 1.0, None)]:
    s = sub[sub.zone == zone]
    axE.scatter(s.astro_z, s.berg_z, s=14 if zone == "border" else 26, c=zone_face[zone],
                edgecolors=ec, linewidths=lw, alpha=0.85, zorder=5 if zone == "border" else 6)
s = ident[ident.type == "Astro-CB"]
axE.scatter(s.astro_z, s.berg_z, s=16, c=C["Astro-CB"], edgecolors="#333333", linewidths=0.4, zorder=5)
s = ident[ident.type == "Bergmann"]
axE.scatter(s.astro_z, s.berg_z, s=30, marker="*", c=C["Bergmann"], edgecolors="#7a1213",
            linewidths=0.3, alpha=0.8, zorder=7)
from matplotlib import patheffects as pe2
for gname, lab, col, off in [("BG_PCL_ref", "BG (PCL) median", GROUP["BG"], (10, 8)),
                             ("parenchymal_ref", "parenchymal astro median", ENTITY["Astro_CB"], (-4, -24))]:
    axE.scatter(isum.loc[gname, "astro_z_med"], isum.loc[gname, "berg_z_med"], s=120, marker="P",
                c=col, edgecolors="#0b0b0b", linewidths=0.8, zorder=9)
    axE.annotate(lab, (isum.loc[gname, "astro_z_med"], isum.loc[gname, "berg_z_med"]),
                 textcoords="offset points", xytext=off, fontsize=7.0, color="#1c1c1c", zorder=10,
                 ha="right" if off[0] < 0 else "left",
                 path_effects=[pe2.withStroke(linewidth=1.8, foreground="white")])
axE.axhline(0, c="#cccccc", lw=0.75, zorder=0); axE.axvline(0, c="#cccccc", lw=0.75, zorder=0)
axE.set_xlabel("astrocyte module (z)"); axE.set_ylabel("Bergmann module (z)")
panel_letter(axE, "e", dx=-0.26, dy=1.06)

# ---------- f: cross-dataset Bergmann share, geometric band -> folial ML ----------
# open bar = raw geometric band; filled bar = folial/ML-proper. The Xenium folial
# denominator is the OUTCOME of testing all 223 band somata, not the swipe pool.
# DO NOT SIMPLIFY (O.G., review 3 point 2): both the two-tier bars with their printed
# denominators (73/223 -> 73/73) and the caption sentence stating that 73/73 is the
# anatomical residue after examining the 223 — not the initial swipe pool — must be
# kept. Losing either turns a defensible result back into a circular-looking one.
W = 0.34
pairs = [(0, Z_RAW, f"{counts['Bergmann']}/{TOT}", Z_PRO, f"{pcounts['Bergmann']}/{PTOT}"),
         (2, 100 * SU_BAND[0] / SU_BAND[1], f"{SU_BAND[0]}/{SU_BAND[1]}", 100.0, f"{n_swipe}/{n_swipe}")]
for x, vr, nr, vf, nf in pairs:
    axF.bar(x - W / 2, vr, W, facecolor="white", edgecolor="#333333", linewidth=0.8)
    axF.bar(x + W / 2, vf, W, color=GROUP["FC"], alpha=0.9, edgecolor="#333333", linewidth=0.6)
    axF.text(x - W / 2, vr + 1.5, nr, ha="center", fontsize=6.9, color="#555555", rotation=90, va="bottom")
    axF.text(x + W / 2, vf + 1.5, nf, ha="center", fontsize=6.9, color="#555555", rotation=90, va="bottom")
axF.bar(1, H_BAND, W, color=GROUP["FC"], alpha=0.9, edgecolor="#333333", linewidth=0.6)
axF.text(1, H_BAND + 1.5, f"{mlbg}/{mlbg + mlnb}", ha="center", fontsize=6.9, color="#555555",
         rotation=90, va="bottom")
from matplotlib.patches import Patch
axF.legend(handles=[Patch(facecolor="white", edgecolor="#333333", label="geometric band"),
                    Patch(facecolor=GROUP["FC"], alpha=0.9, edgecolor="#333333", label="folial / ML-proper")],
           fontsize=6.9, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.46),
           handlelength=1.1, handletextpad=0.4, labelspacing=0.3)
axF.set_xticks(range(3)); axF.set_xticklabels(["MERFISH", "Stereo-seq", "Xenium"], fontsize=7, rotation=90,
                        va="top")
axF.tick_params(axis="x", length=0)
axF.set_ylabel("ML astroglia assigned\nBergmann-class (%)")
axF.set_ylim(0, 132); axF.set_yticks([0, 25, 50, 75, 100])
axF.axhline(100, c="#cccccc", lw=0.75, zorder=0)
panel_letter(axF, "f", dx=-0.30, dy=1.06)

# shared legend row
h = [Line2D([0], [0], marker="*", color="w", markerfacecolor=C["Bergmann"], markeredgecolor="#7a1213",
            markersize=10, label="displaced Bergmann (= FC)"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=C["Astro-NT"], markeredgecolor=NT_EDGE,
            markersize=6, label="Astro-NT in the ML mask (census)"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#c9962e",
            markersize=5, label="Astro-NT outside the census (subclass home: 75% WM)"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=C["Astro-NT"], markeredgecolor="#0b0b0b",
            markersize=7, markeredgewidth=1.1, label="Astro-NT beyond 30 µm (margin/interface)"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=C["Astro-CB"], markeredgecolor="#333333",
            markersize=6, label="parenchymal astrocyte (velate, Astro-CB)"),
     Line2D([0], [0], color=PCL_LINE, lw=1.2, label="PCL line")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=7.0,
           bbox_to_anchor=(0.5, 0.0), columnspacing=0.9, handletextpad=0.4)

save(fig, "Figure2_ML_census")
