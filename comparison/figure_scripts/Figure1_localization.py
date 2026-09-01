#!/usr/bin/env python3
"""Figure 1 — validated FC soma localisation across the three spatial datasets.

The figure is rebuilt directly from the validation/anatomy sources, not from the
expression-analysis tables.  The latter contain deliberate background
subsampling in Hao and must never be used to estimate anatomical prevalence.

Hard failures protect the published quantities:
  * source cell sets are present, exclusive and have the expected cardinality;
  * the Hao numerator is the Mouse 1 manual yes-set and its denominator is the
    complete confident-BG anatomical census across its 32 sections;
  * no valid depth is silently clipped;
  * Hao Mouse 2 is excluded from every active panel and count.
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib import patheffects as pe
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import DS_ORDER, DSL, GROUP, GROUP_LABEL, W2, panel_letter, save, use_style
from project_paths import ROOT as PROJECT_ROOT

use_style()

ROOT = str(PROJECT_ROOT)
TAB = os.path.join(ROOT, "comparison", "tables")
ZH = os.path.join(ROOT, "datasets", "zhang_zhuang_2023_merfish")
HAO = os.path.join(ROOT, "datasets", "hao_2024_stereo_seq")
SU = os.path.join(ROOT, "datasets", "suyama_2026_xenium")

SECTIONS = {
    "zhuang": "Zhuang-ABCA-3.005",
    "hao": "Mouse1_T176",
    "suyama": "Adult",
}
UM_PER_UNIT = {"zhuang": 1000.0, "hao": 10.0, "suyama": 1.0}

LAYER_TINT = {"ML": "#f4f7fb", "PCL": "#f7eadc", "GCL": "#f1eef6", "WM": "#e9edf0"}
PCL_LINE = "#a65a28"
TISSUE_EDGE = "#aeb5bb"
FC_EDGE = "#7a1213"


def _read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


print("== FIGURE 1 SOURCE AUDIT ==")

# ---------------------------------------------------------------------------
# Zhuang: the manually validated FC and canonical BG sets are independent of
# the harmonised expression table.  Distances come from the robust layer model.
z_sets = _read_json(os.path.join(ZH, "validated_analysis", "cell_sets.json"))
z_fc = set(map(str, z_sets["FC_validated"]))
z_bg = set(map(str, z_sets["BG"]))
z_cells = pd.read_parquet(os.path.join(ZH, "cb_cells_multivote.parquet"))
z_cells.index = z_cells.index.astype(str)
z_depth_source = pd.read_parquet(os.path.join(TAB, "robust_layers_zhuang.parquet"))
z_depth_source.index = z_depth_source.index.astype(str)

assert len(z_fc) == 362 and len(z_bg) == 6084   # 4 animals (ABCA-1 integrated 2026-08-18)
assert not (z_fc & z_bg)
assert z_fc <= set(z_cells.index) and z_bg <= set(z_cells.index)
# depth: robust layer model (ABCA-2/3/4) or the ABCA-1 segmentation distance (integrated 2026-08-18)
a1_depth = pd.read_csv(os.path.join(ZH, "layer_segmentation", "BG_assignments_ABCA1.csv"),
                       dtype={"cell_label": str}).set_index("cell_label")["dist_to_PCL_um"]
z_depth_all = z_depth_source["dP_um"].astype(float).combine_first(a1_depth.astype(float))
assert z_fc <= set(z_depth_all.index)
z_depth = z_depth_all.loc[sorted(z_fc)].dropna().to_numpy()
assert len(z_depth) == 362 and np.all(np.isfinite(z_depth))

z_sec = z_cells[z_cells.brain_section_label.eq(SECTIONS["zhuang"])]
z_map_fc = z_sec.loc[z_sec.index.isin(z_fc), ["x", "y"]].copy()
z_map_bg = z_sec.loc[z_sec.index.isin(z_bg), ["x", "y"]].copy()
assert (len(z_map_fc), len(z_map_bg)) == (17, 268)

# ---------------------------------------------------------------------------
# Hao: Mouse 1 is the sole canonical animal. Prevalence uses its full
# segmentation census (7,230 BG), not the BG subsample used for expression.
h_dec_all = pd.read_csv(
    os.path.join(HAO, "swipe_validation", "fc_decisions_hao.csv"),
    dtype={"cell_label": str},
)
h_seg_all = pd.read_csv(
    os.path.join(HAO, "layer_segmentation", "BG_assignments_all.csv"),
    dtype={"cell_label": str},
)
h_dec = h_dec_all[h_dec_all.animal.eq("Mouse1")].copy()
h_seg = h_seg_all[h_seg_all.section.astype(str).str.startswith("Mouse1_")].copy()
h_yes = h_dec[h_dec.decision.str.lower().eq("yes")].copy()
h_bg_rows = h_seg[h_seg.region.eq("BG")].copy()

def _hao_keys(frame):
    return set(zip(frame.section.astype(str), frame.cell_label.astype(str)))

assert h_dec.decision.str.lower().value_counts().to_dict() == {"no": 215, "yes": 89}
assert h_seg.section.nunique() == 32
assert h_seg.region.value_counts().to_dict() == {"BG": 7230, "Ambiguous_BG": 661, "FC": 304}
assert _hao_keys(h_yes) <= _hao_keys(h_seg[h_seg.region.eq("FC")])
assert not (_hao_keys(h_yes) & _hao_keys(h_bg_rows))

h_depth = h_yes.dist_to_PCL_um.astype(float).dropna().to_numpy()
assert len(h_depth) == 89 and np.all(np.isfinite(h_depth))

h_sec_seg = h_seg[h_seg.section.eq(SECTIONS["hao"])]
h_sec_yes = set(h_yes.loc[h_yes.section.eq(SECTIONS["hao"]), "cell_label"])
h_map_fc = h_sec_seg.loc[h_sec_seg.cell_label.isin(h_sec_yes), ["x", "y"]].copy()
h_map_bg = h_sec_seg.loc[h_sec_seg.region.eq("BG"), ["x", "y"]].copy()
assert (len(h_map_fc), len(h_map_bg)) == (12, 118)

# ---------------------------------------------------------------------------
# Suyama: the 44-cell purity + anti-PCL (d_pkj >= 30 um) gate is a subset of the 73 manual yes-calls.  BG
# excludes all 73 yes-calls, preventing rejected-by-purity FC from entering the
# reference denominator.  Distances are the deposited d_pkj geometry values.
s_gate = set(map(str, _read_json(os.path.join(SU, "results", "tables", "suyama_fc_gate.json"))))
s_dec = pd.read_csv(os.path.join(SU, "results", "fc_decisions_suyama.csv"), dtype={"cell_label": str})
s_lam = pd.read_csv(
    os.path.join(SU, "results", "tables", "adult_bergmann_laminar_purkinje.csv"),
    dtype={"cell": str},
)
s_yes = set(s_dec.loc[s_dec.decision.str.lower().eq("yes"), "cell_label"].astype(str))
s_pcl = set(s_lam.loc[s_lam.bgclass.isin(["BG", "BG_PCL"]), "cell"].astype(str))
s_bg = s_pcl - s_yes

assert s_dec.decision.str.lower().value_counts().to_dict() == {"yes": 73, "no": 72}
assert len(s_gate) == 44 and s_gate <= s_yes and s_gate <= set(s_lam.cell)
assert len(s_bg) == 1823 and not (s_gate & s_bg)
s_lam_idx = s_lam.set_index("cell", drop=False)
s_depth = s_lam_idx.loc[sorted(s_gate), "d_pkj"].astype(float).dropna().to_numpy()
assert len(s_depth) == 44 and np.all(np.isfinite(s_depth))
s_map_fc = s_lam[s_lam.cell.isin(s_gate)][["x", "y"]].copy()
s_map_bg = s_lam[s_lam.cell.isin(s_bg)][["x", "y"]].copy()
assert (len(s_map_fc), len(s_map_bg)) == (44, 1823)


FC_N = {"zhuang": len(z_fc), "hao": len(h_yes), "suyama": len(s_gate)}
BG_N = {"zhuang": len(z_bg), "hao": len(h_bg_rows), "suyama": len(s_bg)}
DEPTHS = {"zhuang": z_depth, "hao": h_depth, "suyama": s_depth}
DEPTH_SCOPE = {
    "zhuang": "all validated FC",
    "hao": "all active FC (Mouse 1)",
    "suyama": "all purity-gated FC",
}
MAP_FC = {"zhuang": z_map_fc, "hao": h_map_fc, "suyama": s_map_fc}
MAP_BG = {"zhuang": z_map_bg, "hao": h_map_bg, "suyama": s_map_bg}

EXPECTED = {
    "zhuang": {"share": 362 / 6446 * 100, "median": 91.2},
    "hao": {"share": 89 / 7319 * 100, "median": 86.0},
    "suyama": {"share": 44 / 1867 * 100, "median": 43.3348921380052},
}
for key in DS_ORDER:
    share = 100 * FC_N[key] / (FC_N[key] + BG_N[key])
    median = float(np.median(DEPTHS[key]))
    assert np.isclose(share, EXPECTED[key]["share"], atol=1e-12)
    assert np.isclose(median, EXPECTED[key]["median"], atol=1e-9)
    print(
        f"  {key}: prevalence {FC_N[key]}/{FC_N[key] + BG_N[key]} = {share:.3f}% | "
        f"depth n={len(DEPTHS[key])}, median={median:.1f} µm, "
        f"range={DEPTHS[key].min():.1f}–{DEPTHS[key].max():.1f} µm | "
        f"map={len(MAP_FC[key])} FC/{len(MAP_BG[key])} BG"
    )

# A compact, machine-readable audit table accompanies the figure.
source_rows = []
source_labels = {
    "zhuang": (
        "validated_analysis/cell_sets.json",
        "validated_analysis/cell_sets.json",
        "comparison/tables/robust_layers_zhuang.parquet",
    ),
    "hao": (
        "swipe_validation/fc_decisions_hao.csv (Mouse1)",
        "layer_segmentation/BG_assignments_all.csv (Mouse1)",
        "swipe_validation/fc_decisions_hao.csv (Mouse1)",
    ),
    "suyama": (
        "results/tables/suyama_fc_gate.json",
        "adult_bergmann_laminar_purkinje.csv minus 73 manual yes-calls",
        "results/tables/adult_bergmann_laminar_purkinje.csv",
    ),
}
for key in DS_ORDER:
    fsrc, bsrc, dsrc = source_labels[key]
    d = DEPTHS[key]
    source_rows.append(
        {
            "dataset": key,
            "fc_n": FC_N[key],
            "confident_bg_n": BG_N[key],
            "fc_share_pct": 100 * FC_N[key] / (FC_N[key] + BG_N[key]),
            "depth_n": len(d),
            "depth_scope": DEPTH_SCOPE[key],
            "depth_median_um": np.median(d),
            "depth_min_um": np.min(d),
            "depth_max_um": np.max(d),
            "map_section": SECTIONS[key],
            "map_fc_n": len(MAP_FC[key]),
            "map_bg_n": len(MAP_BG[key]),
            "fc_source": fsrc,
            "bg_source": bsrc,
            "depth_source": dsrc,
        }
    )
pd.DataFrame(source_rows).to_csv(os.path.join(TAB, "Figure1_localization_values.csv"), index=False)


# ---------------------------------------------------------------------------
# Neutral tissue context from the exact display files used during the manual
# swipe.  All context points use one colour: they define the section footprint
# without being mistaken for a second FC/BG or layer classification.
z_swipe = _read_json(os.path.join(ZH, "swipe_validation", "swipe_zhuang.json"))
h_swipe = _read_json(os.path.join(HAO, "swipe_validation", "swipe_data.json"))
s_swipe = pd.read_csv(os.path.join(SU, "results", "swipe_validation", "swipe_cells.csv"))

CONTEXT = {
    "zhuang": pd.DataFrame(
        z_swipe["sections"][SECTIONS["zhuang"]]["cells"],
        columns=["x", "y", "region", "type"],
    )[["x", "y"]],
    "hao": pd.DataFrame(
        h_swipe["sections"][SECTIONS["hao"]]["cells"],
        columns=["x", "y", "region", "type"],
    )[["x", "y"]],
    "suyama": s_swipe[["x", "y"]],
}


def draw_neutral_context(ax, key):
    pts = CONTEXT[key]
    if len(pts) > 7000:
        pts = pts.sample(n=7000, random_state={"zhuang": 101, "hao": 201, "suyama": 301}[key])
    ax.scatter(
        pts.x,
        pts.y,
        s=1.15,
        c="#9aa1a7",
        alpha=0.20,
        edgecolors="none",
        rasterized=True,
        zorder=1,
    )


def map_window(key):
    """Return a common physical half-width and the exact audited-cell midpoint."""
    xy = pd.concat([MAP_FC[key], MAP_BG[key]], ignore_index=True)
    cx = (xy.x.min() + xy.x.max()) / 2
    cy = (xy.y.min() + xy.y.max()) / 2
    half_um = max(xy.x.max() - xy.x.min(), xy.y.max() - xy.y.min()) * UM_PER_UNIT[key] / 2
    return float(half_um), (float(cx), float(cy))


window_needs = {key: map_window(key) for key in DS_ORDER}
COMMON_HALF_UM = max(v[0] for v in window_needs.values()) + 100.0
print(f"  common map half-width: {COMMON_HALF_UM:.0f} µm (identical magnification)")


def draw_schematic(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1)
    ax.axis("off")
    ax.add_patch(Rectangle((0.05, 0.54), 0.90, 0.40, facecolor=LAYER_TINT["ML"], edgecolor="none"))
    ax.add_patch(Rectangle((0.05, 0.46), 0.90, 0.08, facecolor=LAYER_TINT["PCL"], edgecolor="none"))
    ax.add_patch(Rectangle((0.05, 0.08), 0.90, 0.38, facecolor=LAYER_TINT["GCL"], edgecolor="none"))
    ax.plot([0.05, 0.95], [0.94, 0.94], color="#7f8791", lw=0.8)
    ax.plot([0.05, 0.95], [0.50, 0.50], color=PCL_LINE, lw=1.2)

    # Both cells share a radial Bergmann morphology; only soma position differs.
    ax.plot([0.37, 0.37], [0.51, 0.92], color=GROUP["BG"], lw=1.0, alpha=0.70)
    ax.plot([0.68, 0.68], [0.73, 0.92], color=GROUP["FC"], lw=1.0, alpha=0.70)
    ax.scatter([0.37], [0.50], s=40, c=GROUP["BG"], edgecolors="#26323b", linewidths=0.5, zorder=4)
    ax.scatter([0.68], [0.72], s=100, marker="*", c="white", edgecolors="none", zorder=4)
    ax.scatter([0.68], [0.72], s=52, marker="*", c=GROUP["FC"], edgecolors=FC_EDGE, linewidths=0.5, zorder=5)

    halo = [pe.withStroke(linewidth=2.0, foreground="white")]
    for label, y in [("ML", 0.75), ("PCL", 0.50), ("GCL", 0.27)]:
        ax.text(0.10, y, label, ha="left", va="center", fontsize=7, color="#5a5a57", path_effects=halo)
    ax.text(0.37, 0.41, "BG soma", ha="center", va="top", fontsize=7.5, color=GROUP["BG"])
    ax.text(0.68, 0.62, "FC soma", ha="center", va="top", fontsize=7.5, color=GROUP["FC"])
    ax.text(0.50, 0.98, "pia", ha="center", va="bottom", fontsize=7, color="#5a5a57")


def add_scale_bar(ax, key):
    x0, x1 = ax.get_xlim()
    y0, y1 = sorted(ax.get_ylim())
    length = 500.0 / UM_PER_UNIT[key]
    x_start = x1 - 0.07 * (x1 - x0) - length
    y = y0 + 0.055 * (y1 - y0)
    ax.plot([x_start, x_start + length], [y, y], color="#111111", lw=1.6, solid_capstyle="butt", zorder=8)
    ax.text(
        x_start + length / 2,
        y + 0.025 * (y1 - y0),
        "500 µm",
        ha="center",
        va="bottom",
        fontsize=7,
        path_effects=[pe.withStroke(linewidth=2.0, foreground="white")],
        zorder=9,
    )


# ---------------------------------------------------------------------------
# Figure layout: definition + three same-scale maps; exact prevalence and all
# valid depths below.  No inferential CI is attached to cell counts.
fig = plt.figure(figsize=(W2, 5.25))
gs = GridSpec(
    3,
    16,                     # 16 cols: one empty gutter column between e|f and between f|g
    figure=fig,
    height_ratios=[0.54, 0.08, 0.46],
    hspace=0.20,
    wspace=0.85,
    left=0.07,
    right=0.985,
    bottom=0.11,
    top=0.93,
)
ax_schematic = fig.add_subplot(gs[0, 0:3])
ax_maps = [fig.add_subplot(gs[0, 4:8]), fig.add_subplot(gs[0, 8:12]), fig.add_subplot(gs[0, 12:16])]
ax_legend = fig.add_subplot(gs[1, :])
ax_share = fig.add_subplot(gs[2, 0:6])
ax_depth = fig.add_subplot(gs[2, 7:11])
ax_norm = fig.add_subplot(gs[2, 12:16])

draw_schematic(ax_schematic)
panel_letter(ax_schematic, "a", dx=0.00, dy=1.10)

map_labels = {
    "zhuang": "ABCA-3.005",
    "hao": "Mouse1_T176",
    "suyama": "adult section",
}
# Datasets are labelled by PLATFORM throughout the figures (O.G., review 3 point 2); the first
# figure carries the first author once, in parentheses, and every later figure is platform-only.
# Note "Zhang" and not "Zhuang": the MERFISH atlas is cited as Zhang et al. 2023 ("Zhuang" is the
# laboratory and survives only inside data identifiers such as Zhuang-ABCA-3.005).
DSL_FIG1 = {"zhuang": "MERFISH (Zhang)", "hao": "Stereo-seq (Hao)", "suyama": "Xenium (Suyama)"}
for i, key in enumerate(DS_ORDER):
    ax = ax_maps[i]
    draw_neutral_context(ax, key)

    bg = MAP_BG[key]
    fc = MAP_FC[key]
    bg_size = {"zhuang": 4.3, "hao": 5.2, "suyama": 2.0}[key]
    ax.scatter(bg.x, bg.y, s=bg_size, c=GROUP["BG"], alpha=0.58, edgecolors="none", zorder=6)
    ax.scatter(fc.x, fc.y, s=95, marker="*", c="white", edgecolors="none", zorder=6)
    ax.scatter(
        fc.x,
        fc.y,
        s=48,
        marker="*",
        c=GROUP["FC"],
        edgecolors=FC_EDGE,
        linewidths=0.4,
        zorder=7,
    )

    _, (cx, cy) = window_needs[key]
    half = COMMON_HALF_UM / UM_PER_UNIT[key]
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal")
    ax.set_box_aspect(1)
    ax.axis("off")
    ax.set_title(
        f"{DSL_FIG1[key]}\n{map_labels[key]} · {len(fc)} FC / {len(bg)} BG",
        fontsize=8,
        fontweight="normal",
        linespacing=1.35,
        pad=4,
    )
    add_scale_bar(ax, key)
    panel_letter(ax, "bcd"[i], dx=-0.04, dy=1.18)

ax_legend.axis("off")
legend_handles = [
    Line2D(
        [0], [0], marker="*", linestyle="none", markerfacecolor=GROUP["FC"],
        markeredgecolor=FC_EDGE, markersize=10, label="validated FC"
    ),
    Line2D(
        [0], [0], marker="o", linestyle="none", markerfacecolor=GROUP["BG"],
        markeredgecolor="none", alpha=0.70, markersize=5.5, label="confident BG"
    ),
]
ax_legend.legend(
    handles=legend_handles,
    loc="center",
    ncol=2,
    frameon=False,
    fontsize=7.5,
    handletextpad=0.45,
    columnspacing=1.5,
)

# Panel e — exact observed proportions.  No cell-level confidence interval:
# cells from the same section/animal are not independent biological replicates.
ypos = np.arange(3)
shares = np.array([100 * FC_N[k] / (FC_N[k] + BG_N[k]) for k in DS_ORDER])
for y, value in zip(ypos, shares):
    ax_share.plot([0, value], [y, y], color="#d6d8da", lw=2.2, solid_capstyle="round", zorder=1)
ax_share.scatter(shares, ypos, s=48, c=GROUP["FC"], edgecolors=FC_EDGE, linewidths=0.6, zorder=3)
# Biological replicates as individual points (review point 4, 2026-08-17):
# small open circles = MERFISH animals (row 0) and Stereo-seq FC-containing sections (row 1).
zc_rep = pd.read_parquet(os.path.join(TAB, "zhuang_astro_cells.parquet"))
per_an = zc_rep.groupby("animal")["group"].value_counts().unstack().fillna(0)
an_shares = np.sort((100 * per_an["FC"] / (per_an["FC"] + per_an["BG"])).values)
assert np.allclose(np.round(an_shares, 2), [3.25, 3.85, 5.23, 6.40], atol=0.01), an_shares
sec_tab = pd.read_csv(os.path.join(TAB, "hao_fc_share_by_section.csv"))
sec_shares = sec_tab.loc[sec_tab.FC > 0, "share"].values
assert len(sec_shares) == 26 and abs(np.median(sec_shares) - 1.45) < 0.05, (len(sec_shares), np.median(sec_shares))
rng_rep = np.random.default_rng(4)
ax_share.scatter(an_shares, 0 + rng_rep.uniform(-0.14, 0.14, len(an_shares)),
                 s=16, facecolors="white", edgecolors=FC_EDGE, linewidths=0.8, zorder=4)
ax_share.scatter(sec_shares, 1 + rng_rep.uniform(-0.16, 0.16, len(sec_shares)),
                 s=9, facecolors="white", edgecolors=FC_EDGE, linewidths=0.55, alpha=0.8, zorder=4)
ax_share.text(0.99, 0.97, "open circles: individual animals (MERFISH) / sections (Stereo-seq)",
              transform=ax_share.transAxes, ha="right", va="top", fontsize=7.0, color="#606060")
for y, key, value in zip(ypos, DS_ORDER, shares):
    denominator = FC_N[key] + BG_N[key]
    ax_share.text(value + 0.15, y - 0.16, f"{value:.2f}%", ha="left", va="center", fontsize=8,
                  fontweight="bold", zorder=7,
                  path_effects=[pe.withStroke(linewidth=2.2, foreground="white")])
    ax_share.text(
        value + 0.15,
        y + 0.32,
        f"{FC_N[key]} / {denominator}",
        ha="left",
        va="center",
        fontsize=7,
        color="#606060",
    )
ax_share.set_yticks(ypos)
ax_share.set_yticklabels(["MERFISH", "Stereo-seq", "Xenium"], fontsize=8)
ax_share.invert_yaxis()
ax_share.set_xlim(0, 9.8)
ax_share.set_xticks([0, 2, 4, 6, 8])
ax_share.set_xlabel("FC / (FC + confident BG)  (%)")
ax_share.tick_params(axis="y", length=0)
ax_share.spines[["top", "right", "left"]].set_visible(False)
panel_letter(ax_share, "e", dx=-0.17, dy=1.08)

# Panel f — every physically calibrated depth, without the old <=400-µm clip.
positions = np.arange(3)
depth_values = [DEPTHS[k] for k in DS_ORDER]
ax_depth.boxplot(
    depth_values,
    positions=positions,
    widths=0.48,
    patch_artist=True,
    showfliers=False,
    showcaps=True,
    medianprops={"color": "#111111", "lw": 1.5},
    boxprops={"facecolor": "#f8dede", "edgecolor": "#4b4b4b", "lw": 0.9},
    whiskerprops={"color": "#4b4b4b", "lw": 0.9},
    capprops={"color": "#4b4b4b", "lw": 0.9},
)
rng = np.random.default_rng(20260812)
for x, values in zip(positions, depth_values):
    jitter = rng.uniform(-0.12, 0.12, len(values))
    ax_depth.scatter(
        x + jitter,
        values,
        s=7,
        c=GROUP["FC"],
        alpha=0.46,
        edgecolors="none",
        rasterized=False,
        zorder=3,
    )
# per-dataset median values are stated in the legend (panel too narrow for in-panel labels)
ax_depth.set_xticks(positions)
ax_depth.set_xticklabels([f"MERFISH\nn={len(DEPTHS['zhuang'])}", f"Stereo-seq\nM1, n={len(DEPTHS['hao'])}",
                          f"Xenium\nn={len(DEPTHS['suyama'])}"], fontsize=7, rotation=18, ha="right",
                         rotation_mode="anchor")
ax_depth.tick_params(axis="x", length=0)
ax_depth.set_ylabel("soma distance to PCL (µm)")
ax_depth.set_ylim(0, 800)
ax_depth.set_yticks([0, 200, 400, 600, 800])
ax_depth.spines[["top", "right"]].set_visible(False)
panel_letter(ax_depth, "f", dx=-0.16, dy=1.20)

# Panel g — normalized ML depth, 0 = PCL, 1 = pia (review point 4, 2026-08-17).
# Zhuang: exact per-soma depth/(depth+pia); Hao & Suyama: depth / local ML thickness
# (95th pct within 100 um; build_normalized_depth.py).
ndep = pd.read_csv(os.path.join(TAB, "fc_normalized_depth.csv"))
nd_vals = [ndep.loc[ndep.dataset == k, "norm_depth"].values for k in DS_ORDER]
nd_med = [float(np.median(v)) for v in nd_vals]
assert np.allclose(nd_med, [0.572, 0.527, 0.230], atol=0.01), nd_med
ax_norm.boxplot(nd_vals, positions=positions, widths=0.48, patch_artist=True, showfliers=False,
                showcaps=True, medianprops={"color": "#111111", "lw": 1.5},
                boxprops={"facecolor": "#f8dede", "edgecolor": "#4b4b4b", "lw": 0.9},
                whiskerprops={"color": "#4b4b4b", "lw": 0.9}, capprops={"color": "#4b4b4b", "lw": 0.9})
rng_n = np.random.default_rng(20260817)
for xpos, values in zip(positions, nd_vals):
    ax_norm.scatter(xpos + rng_n.uniform(-0.12, 0.12, len(values)), values, s=7,
                    c=GROUP["FC"], alpha=0.46, edgecolors="none", zorder=3)
ax_norm.axhline(0, color="#a65a28", lw=0.9)
ax_norm.axhline(1, color="#898781", lw=0.9, ls=(0, (4, 3)))
ax_norm.text(2.62, 0.015, "PCL", fontsize=7, color="#a65a28", ha="right", va="bottom")
ax_norm.text(2.62, 0.985, "pia", fontsize=7, color="#898781", ha="right", va="top")
ax_norm.set_xticks(positions)
ax_norm.set_xticklabels([f"MERFISH\nn={len(nd_vals[0])}", f"Stereo-seq\nM1, n={len(nd_vals[1])}",
                         f"Xenium\nn={len(nd_vals[2])}"], fontsize=7, rotation=18, ha="right",
                        rotation_mode="anchor")
ax_norm.tick_params(axis="x", length=0)
ax_norm.set_ylabel("normalized ML depth\n(0 = PCL, 1 = pia)")
ax_norm.set_ylim(-0.05, 1.10)
ax_norm.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_norm.spines[["top", "right"]].set_visible(False)
panel_letter(ax_norm, "g", dx=-0.24, dy=1.20)

save(fig, "Figure1_localization")
