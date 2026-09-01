#!/usr/bin/env python
"""Figure S13 — FC share by cerebellar lobule / territory (review point 18).

a: FC/(FC + confident BG) per Allen-CCF structure, anatomical order (vermal
   lobules I-X, then hemispheric SIM/AN/PRM/COPY, then FL/PFL); n annotated.
b: aggregate territories. The profile is a within-animal observation: ABCA-3
   samples all 14 structures; the two sparse animals contribute counts at the
   mediolateral levels they cover (caption).
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import use_style, GROUP, W2, panel_letter, save

use_style()
from project_paths import ROOT as PROJECT_ROOT
TAB = str(PROJECT_ROOT) + "/comparison/tables"

ORDER = ["LING", "CENT", "CUL", "DEC", "FOTU", "PYR", "UVU", "NOD",
         "SIM", "AN", "PRM", "COPY", "FL", "PFL"]
ROMAN = {"LING": "I", "CENT": "II-III", "CUL": "IV-V", "DEC": "VI", "FOTU": "VII",
         "PYR": "VIII", "UVU": "IX", "NOD": "X"}
TERR_ORDER = ["anterior vermis (I-V)", "posterior vermis (VI-IX)", "nodulus (X)",
              "hemispheres", "flocculus/paraflocculus"]

t = pd.read_csv(TAB + "/fc_share_by_lobule.csv")
st = t[(t.level == "structure") & (t.animal == "all")].set_index("unit").loc[ORDER]
assert int(st.n_FC.sum()) == 351 and int(st.n_BG.sum()) == 6055  # 4 animals; 11 FC / 29 BG CB-unassigned
assert abs(st.loc["NOD", "share_pct"] - 14.80) < 0.01
assert abs(st.loc["UVU", "share_pct"] - 7.73) < 0.01

fig, (axA, axB) = plt.subplots(1, 2, figsize=(W2, 3.6),
                               gridspec_kw=dict(width_ratios=[1.15, 1], wspace=0.85,
                                                top=0.87, bottom=0.17))

# ---- a: per structure (horizontal bars, anatomical order top -> bottom) ----
y = np.arange(len(ORDER))[::-1]
axA.barh(y, st.share_pct, height=0.62, color=GROUP["FC"], edgecolor="#333333", lw=0.4)
for yi, (unit, r) in zip(y, st.iterrows()):
    axA.text(r.share_pct + 0.35, yi, f"{int(r.n_FC)}", va="center", fontsize=7,
             color="#52514e")
axA.axhline(y[7] - 0.5, color="#c3c2b7", lw=0.9)
axA.axhline(y[11] - 0.5, color="#c3c2b7", lw=0.9)
axA.text(21.8, y[3] - 0.4, "vermis\n(I-X)", ha="right", fontsize=7.5, color="#52514e")
axA.text(21.8, y[9] - 0.6, "hemispheres", ha="right", fontsize=7.5, color="#52514e")
axA.text(21.8, y[12] - 0.7, "flocculus/\nparaflocculus", ha="right", fontsize=7.5,
         color="#52514e")
axA.set_yticks(y)
axA.set_yticklabels([f"{u} {ROMAN[u]}" if u in ROMAN else u for u in ORDER], fontsize=7.5)
axA.set_xlabel("FC share of Bergmann-class cells (%)\nFC / (FC + confident BG); n FC at bar end")
axA.set_xlim(0, 22.5)
axA.set_ylim(-0.7, len(ORDER) - 0.3)
axA.set_title("FC share by Allen-CCF structure", loc="left",
              fontweight="bold", fontsize=9.5)
panel_letter(axA, "a", dx=-0.20, dy=1.10)

# ---- b: aggregate territories ----
te = t[(t.level == "territory") & (t.animal == "all")].set_index("unit").loc[TERR_ORDER]
yb = np.arange(len(TERR_ORDER))[::-1]
axB.barh(yb, te.share_pct, height=0.6, color=GROUP["FC"], edgecolor="#333333", lw=0.4)
for yi, (unit, r) in zip(yb, te.iterrows()):
    axB.text(r.share_pct + 0.35, yi, f"{int(r.n_FC)}/{int(r.n_FC + r.n_BG)}",
             va="center", fontsize=7, color="#52514e")
axB.set_yticks(yb)
axB.set_yticklabels(["anterior vermis (I-V)", "posterior vermis (VI-IX)", "nodulus (X)",
                     "hemispheres", "flocculus/parafloccul."], fontsize=7.5)
axB.set_xlabel("FC share (%)")
axB.set_xlim(0, 23.5)
axB.set_title("Territories", loc="left", fontweight="bold", fontsize=9.5)
panel_letter(axB, "b", dx=-0.62, dy=1.10)

save(fig, "FigureS13_lobules")
