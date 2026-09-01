#!/usr/bin/env python
"""Figure 2 — molecular-layer astroglial census in Zhuang MERFISH.

Panels show class-first composition, white-matter distance, and the spatial
census map on section ABCA-3.005. Hao controls are reported in Figure S2.
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

MAP_SEC = "Zhuang-ABCA-3.005"   # same section as Figure 1; census coverage 100% (G5 sweep 2026-08-10)
C = {"Bergmann": GROUP["FC"], "Astro-NT": ENTITY["Astro_NT"], "Astro-CB": ENTITY["Astro_CB"]}
NT_PT, NT_EDGE = "#d08c00", "#7a5200"
LAYER_TINT = {"ML": "#dcebf9", "PCL": "#f6d9a8", "GCL": "#e7ddf3", "WM": "#e4e7eb"}
PCL_LINE = "#b8560f"
ORDER3 = ["Bergmann", "Astro-NT", "Astro-CB"]


def wilson(k, n, z=1.959964):
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    hw = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return 100 * (ctr - hw), 100 * (ctr + hw)


# ================= G1 — DATA AUDIT (aborts on mismatch) =================
print("== DATA AUDIT (Figure 2) ==")
zc = pd.read_parquet(TAB + "/zhuang_ML_census_cells.parquet")
z45 = zc[zc.depth_from_PCL >= 45]
counts = {t: int((z45.type == t).sum()) for t in ORDER3}
assert counts == {"Bergmann": 263, "Astro-NT": 246, "Astro-CB": 10}, counts
TOT = sum(counts.values()); assert TOT == 519
meds = {t: float(z45.loc[z45.type == t, "dist_to_WM"].median()) for t in ORDER3}
assert round(meds["Bergmann"]) == 96 and round(meds["Astro-NT"]) == 18 and round(meds["Astro-CB"]) == 66, meds
frac30 = float((z45.loc[z45.type == "Astro-NT", "dist_to_WM"] <= 30).mean())
assert round(100 * frac30) == 72, frac30
pia_med = float(z45.loc[z45.type == "Bergmann", "dist_to_pia"].median())
assert round(pia_med) == 40, pia_med
lbt = pd.read_csv(TAB + "/zhuang_layer_by_type.csv", index_col=0)
nt = lbt.loc[[i for i in lbt.index if "Astro-NT" in str(i) or "318" in str(i)][0]]
cb = lbt.loc[[i for i in lbt.index if "Astro-CB" in str(i) or "317" in str(i)][0]]
ntw = nt["WM"] / nt.sum(); cbg = cb["GCL"] / cb.sum(); cbm = cb["ML"] / cb.sum()
assert round(100 * ntw) == 75 and round(100 * cbg) == 85 and round(100 * cbm, 1) == 1.2
gmap = z45[z45.section == MAP_SEC]
mcounts = {t: int((gmap.type == t).sum()) for t in ORDER3}
assert mcounts == {"Bergmann": 18, "Astro-NT": 15, "Astro-CB": 3}, mcounts
print(f"  census: 263/246/10 (total {TOT}) ✓  WM med 96/18/66 ✓  72%<30µm ✓  pia(displaced)=40µm ✓")
print(f"  upstream: Astro-NT {100*ntw:.1f}% in WM ✓  Astro-CB {100*cbg:.1f}% GCL / {100*cbm:.2f}% ML ✓")
print(f"  map section {MAP_SEC}: {mcounts} ✓")

# map region + coverage check
zz = pd.read_parquet(ZH + "/cb_cells_multivote.parquet")
zz = zz[zz.brain_section_label == MAP_SEC]
reg = layer_regions(zz.x, zz.y, zz.is_gran, zz.subclass.astype(str).str.contains("Oligo"), umpp=1000.0)
XX, YY, T = reg["XX"], reg["YY"], np.asarray(reg["tissue"])
tree = cKDTree(np.column_stack([XX.ravel(), YY.ravel()]))
_, flat = tree.query(np.column_stack([gmap.x.values, gmap.y.values]))
cov = float(T.ravel()[flat].mean())
assert cov == 1.0, f"census coverage on {MAP_SEC} = {cov}"
print("  map coverage of census cells: 100% ✓")

# ================= figure =================
fig = plt.figure(figsize=(W2, 2.95))
gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.15, 1.12], wspace=0.56,
              top=0.90, bottom=0.24, left=0.08, right=0.985)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[0, 2])

# ---------- a: composition (%) with Wilson CIs ----------
x = np.arange(3)
pct = [100 * counts[t] / TOT for t in ORDER3]
ci = [wilson(counts[t], TOT) for t in ORDER3]
err = np.array([[v - lo, hi - v] for v, (lo, hi) in zip(pct, ci)]).T
axA.bar(x, pct, width=0.52, color=[C[t] for t in ORDER3], alpha=0.9,
        edgecolor="#333333", linewidth=0.6)
axA.errorbar(x, pct, yerr=err, fmt="none", ecolor="#0b0b0b", elinewidth=0.9, capsize=2.5)
axA.set_xticks(x); axA.set_xticklabels(["displaced\nBergmann", "WM\nastro", "parench.\nastro"], fontsize=7.5)
axA.tick_params(axis="x", length=0)
axA.set_ylabel("share of ML-mask\nastroglia (%)")
axA.set_ylim(0, 62); axA.set_yticks([0, 20, 40, 60])
panel_letter(axA, "a", dx=-0.34, dy=1.10)

# ---------- b: soma-to-WM distance ----------
rng = np.random.default_rng(0)
dmax = 0
for i, t in enumerate(ORDER3):
    d = z45.loc[z45.type == t, "dist_to_WM"].values
    dmax = max(dmax, d.max())
    if len(d) >= 25:
        pc = axB.violinplot([d], positions=[i], showextrema=False, widths=0.75)["bodies"][0]
        pc.set_facecolor(C[t]); pc.set_alpha(0.15)
    fill = NT_PT if t == "Astro-NT" else C[t]
    edge = NT_EDGE if t == "Astro-NT" else "none"
    axB.scatter(i + rng.normal(0, 0.05, len(d)), d, s=5, c=fill, alpha=0.5,
                edgecolors=edge, linewidths=0.3)
    axB.hlines(np.median(d), i - 0.3, i + 0.3, color="#0b0b0b", lw=1.3, zorder=5)
axB.axhline(30, ls=(0, (4, 3)), c="#898781", lw=0.9)
axB.set_xticks(range(3)); axB.set_xticklabels(["displaced\nBergmann", "WM\nastro", "parench.\nastro"], fontsize=7.5)
axB.tick_params(axis="x", length=0)
axB.set_xlim(-0.6, 2.6); axB.set_ylim(0, min(dmax * 1.05, 430))
axB.set_ylabel("soma-to-white-matter\ndistance (µm)")
axB.set_yticks([0, 100, 200, 300, 400])
panel_letter(axB, "b", dx=-0.30, dy=1.10)

# ---------- c: spatial census map (section ABCA-3.005, as in Figure 1) ----------
paint_layers(axC, reg, colors=LAYER_TINT, pcl_line=PCL_LINE, pcl_lw=0.9)
Ts = ndi.gaussian_filter(T.astype(float), sigma=1.6)
axC.contour(XX, YY, Ts, levels=[0.5], colors=["#b6bdc4"], linewidths=0.55, zorder=8)
for t in ["Astro-NT", "Astro-CB"]:
    gg = gmap[gmap.type == t]
    ec = "#7a5200" if t == "Astro-NT" else "#333333"
    axC.scatter(gg.x, gg.y, s=16, c=C[t], edgecolors=ec, linewidths=0.45, zorder=6)
gg = gmap[gmap.type == "Bergmann"]   # displaced Bergmann = the FC star of Figure 1 (same glyph)
axC.scatter(gg.x, gg.y, s=80, marker="*", c="white", edgecolors="none", zorder=5)
axC.scatter(gg.x, gg.y, s=40, marker="*", c=C["Bergmann"], edgecolors="#7a1213",
            linewidths=0.45, alpha=0.95, zorder=6)
pad = 0.15
x0, x1 = gmap.x.min() - pad, gmap.x.max() + pad
y0, y1 = gmap.y.min() - pad, gmap.y.max() + pad
cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
half = max(x1 - x0, y1 - y0) / 2
axC.set_xlim(cx - half, cx + half); axC.set_ylim(cy - half, cy + half)
axC.set_aspect("equal"); axC.set_box_aspect(1); axC.axis("off")
xl0, xl1 = axC.get_xlim(); yl0, yl1 = sorted(axC.get_ylim())
L = 0.5
xs = xl1 - (xl1 - xl0) * 0.05 - L; ys = yl0 + (yl1 - yl0) * 0.05
axC.add_line(Line2D([xs, xs + L], [ys, ys], color="#0b0b0b", lw=1.8, solid_capstyle="butt", zorder=11))
from matplotlib import patheffects as pe
axC.text(xs + L / 2, ys + (yl1 - yl0) * 0.022, "500 µm", ha="center", va="bottom", fontsize=7,
         zorder=11, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
# direct layer labels, window- and cell-aware (same rule as Figure 1)
sd = reg["sdist"]; wmm = reg["wmask"] & T
mx, my = 0.07 * (xl1 - xl0), 0.07 * (yl1 - yl0)
inwin = (XX > xl0 + mx) & (XX < xl1 - mx) & (YY > yl0 + my) & (YY < yl1 - my)
rmin = 0.06 * (xl1 - xl0)
far = np.ones_like(XX, dtype=bool)
for sx_, sy_ in zip(gmap.x.values, gmap.y.values):
    far &= ((XX - sx_) ** 2 + (YY - sy_) ** 2) > rmin ** 2
ok = inwin & far
spots = {}
mlf = np.where(T & ~wmm & ok, sd, np.nan)
if np.isfinite(mlf).any():
    spots["ML"] = np.unravel_index(np.nanargmin(mlf), mlf.shape)
    spots["GCL"] = np.unravel_index(np.nanargmax(mlf), mlf.shape)
if (wmm & ok).sum() > 20:
    dt = ndi.distance_transform_edt(wmm & ok)
    spots["WM"] = np.unravel_index(np.argmax(dt), dt.shape)
for name, (ii2, jj2) in spots.items():
    axC.text(XX[ii2, jj2], YY[ii2, jj2], name, fontsize=7, color="#5a5a57", ha="center",
             va="center", zorder=9, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
panel_letter(axC, "c", dx=0.0, dy=1.10)

# shared legend row (serves all three panels)
h = [Line2D([0], [0], marker="*", color="w", markerfacecolor=C["Bergmann"],
            markeredgecolor="#7a1213", markersize=10, label="displaced Bergmann (= FC)"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=C["Astro-NT"],
            markeredgecolor="#7a5200", markersize=6, label="WM astrocyte (Astro-NT)"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=C["Astro-CB"],
            markeredgecolor="#333333", markersize=6, label="parenchymal astrocyte (velate, Astro-CB)"),
     Line2D([0], [0], color=PCL_LINE, lw=1.2, label="PCL line")]
fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, fontsize=7.5,
           bbox_to_anchor=(0.5, 0.0), columnspacing=1.3, handletextpad=0.45)

save(fig, "Figure2_ML_census")
