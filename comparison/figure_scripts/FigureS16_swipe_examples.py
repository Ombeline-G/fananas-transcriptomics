#!/usr/bin/env python
"""Figure S16 — representative accepted and rejected swipe candidates
(review point 26).

Panels are re-rendered from the exact data shown to the rater
(swipe_zhuang.json: per-section cell maps colored by segmentation layer, BG
and Purkinje somata highlighted, candidate soma as the red star; 250-um zoom,
as in the tool). The rater saw ONLY this geometry — no expression values —
so the review is blinded to transcriptomic profile by construction (Methods).
Selection of examples is deterministic: seed-0 sample of 4 accepted and 4
rejected candidates from fc_decisions_zhuang.csv, with two named substitutions
in the rejected row (see SWAPS below).
"""
import json
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT, ZHUANG

SWIPE = ZHUANG / "swipe_validation"
LAYER_DOT = {0: "#9cc6e6", 1: "#c4a3da", 2: "#aab2bb", 3: "#e0b483"}  # ML/GCL/WM/PCL band
C_PURK, C_BGPCL, C_STAR = "#e8894a", "#15348f", "#d62728"
ZOOM_MM = 0.250

data = json.load(open(SWIPE / "swipe_zhuang.json"))
dec = pd.read_csv(SWIPE / "fc_decisions_zhuang.csv", dtype={"cell_label": str})
assert dec.decision.value_counts().to_dict() == {"yes": 190, "no": 50}
cards = pd.DataFrame(data["cards"])
cards["id"] = cards.id.astype(str)
cards = cards.merge(dec[["cell_label", "decision"]], left_on="id", right_on="cell_label")
assert len(cards) == 240

rng = np.random.default_rng(0)
# Accepted: seed-0 sample. Rejected: seed-0 sample with TWO named substitutions, at
# the rater's request (O.G., review 2 point 13; these supersede the single ABCA-3.004
# substitution of review 1): the second and third draws (ABCA-4.003 at 80 um and
# ABCA-4.003 at 97 um) are replaced by candidates whose neighbourhoods are densely
# sampled, so that the geometry behind the call is visible — ABCA-3.017 (49 um above
# the PCL; the densest rejected window, 430 cells; nearest Purkinje-layer soma 54 um)
# and ABCA-4.003 (49 um; 21 Purkinje-layer somata within 200 um, nearest 49 um;
# measured, not asserted). This is a legibility choice, NOT a demonstration of a
# threshold: the manual call is a judgment on the whole neighbourhood, and an ACCEPTED
# panel in this same figure (ABCA-3.015) has its nearest Purkinje-layer soma at 45 um.
# The substitutions are explicit and every one of the 240 decisions is released.
SWAPS = {"49615121079511868065563556194399505935":   # draw 2 -> ABCA-3.017
         "208353139449262393390543151337108120440",
         "247287723180358287728005649534653069466":   # draw 3 -> ABCA-4.003
         "322556046161549581757114812417145617866"}
no_pool = cards[cards.decision == "no"]
no_pick = no_pool.sample(4, random_state=0)
assert set(SWAPS) < set(no_pick.id) and set(SWAPS.values()) < set(no_pool.id)
assert not set(SWAPS.values()) & set(no_pick.id)
no_pick = pd.concat([cards[cards.id == SWAPS[cid]] if cid in SWAPS
                     else no_pick[no_pick.id == cid] for cid in no_pick.id])
picks = {"yes": cards[cards.decision == "yes"].sample(4, random_state=0), "no": no_pick}

fig, axes = plt.subplots(2, 4, figsize=(W2, 4.4),
                         gridspec_kw=dict(wspace=0.12, hspace=0.30, top=0.90, bottom=0.03))
for row, (key, label) in enumerate([("yes", "accepted"), ("no", "rejected")]):
    for col, (_, card) in enumerate(picks[key].iterrows()):
        ax = axes[row, col]
        sec = data["sections"][card.section]
        cells = np.array(sec["cells"], dtype=float)
        cx, cy = card.cx, card.cy
        win = ((np.abs(cells[:, 0] - cx) < ZOOM_MM) & (np.abs(cells[:, 1] - cy) < ZOOM_MM))
        w = cells[win]
        for lay, color in LAYER_DOT.items():
            mm = (w[:, 2] == lay) & (w[:, 3] == 0)
            ax.scatter(w[mm, 0], w[mm, 1], s=5, c=color, linewidths=0)
        # role 1 = Purkinje, role 2 = BG (gen_swipe_from_segmentation.py); the swipe tool
        # draws Purkinje orange and BG_PCL blue — this figure claims to reproduce it exactly
        for role, color, size in [(2, C_BGPCL, 8), (1, C_PURK, 12)]:
            mm = w[:, 3] == role
            ax.scatter(w[mm, 0], w[mm, 1], s=size, c=color, linewidths=0)
        ax.scatter([cx], [cy], marker="*", s=170, c=C_STAR, edgecolors="#000000",
                   linewidths=0.8, zorder=6)
        ax.set_xlim(cx - ZOOM_MM, cx + ZOOM_MM)
        ax.set_ylim(cy - ZOOM_MM, cy + ZOOM_MM)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect("equal")
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_color("#2ca02c" if key == "yes" else "#d62728")
            sp.set_linewidth(1.4)
        ax.set_title(f"{card.section.replace('Zhuang-', '')} · {int(card.dist_pcl)} µm",
                     fontsize=7.5, pad=2)
    axes[row, 0].set_ylabel(label, fontsize=10, fontweight="bold")
from matplotlib.lines import Line2D
key = [Line2D([], [], ls="none", marker="o", ms=4, mfc=LAYER_DOT[0], mec="none", label="ML"),
       Line2D([], [], ls="none", marker="o", ms=4, mfc=LAYER_DOT[3], mec="none", label="PCL band"),
       Line2D([], [], ls="none", marker="o", ms=4, mfc=LAYER_DOT[1], mec="none", label="GCL"),
       Line2D([], [], ls="none", marker="o", ms=4, mfc=LAYER_DOT[2], mec="none", label="WM"),
       Line2D([], [], ls="none", marker="o", ms=5, mfc=C_PURK, mec="none", label="Purkinje soma"),
       Line2D([], [], ls="none", marker="o", ms=4, mfc=C_BGPCL, mec="none", label="BG soma"),
       Line2D([], [], ls="none", marker="*", ms=9, mfc=C_STAR, mec="#000000", mew=0.6, label="candidate")]
fig.legend(handles=key, loc="lower center", ncol=7, frameon=False, fontsize=7.5,
           handletextpad=0.3, columnspacing=1.1, bbox_to_anchor=(0.5, -0.005))
panel_letter(axes[0, 0], "a", dx=-0.10, dy=1.22)
fig.suptitle("Swipe review, exactly as seen by the rater — geometry only, no expression "
             "(MERFISH; 500-µm windows)", fontsize=9.5, fontweight="bold", y=0.985)
save(fig, "FigureS16_swipe_examples")
