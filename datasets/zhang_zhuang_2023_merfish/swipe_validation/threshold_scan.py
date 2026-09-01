#!/usr/bin/env python
"""
Objectify the N <-> purity trade-off of the ML_MIN_DEPTH (45 um) gate.

Two complementary views, both anchored on dist_to_PCL_um:

(A) PURITY vs distance, using the swipe GOLD labels.
    The 240 FC candidates were all called by segmentation at dist>=45 um, then
    curated by hand (yes/no). For a candidate threshold D, among candidates with
    dist>=D we report N (cells kept), N_yes (human-confirmed), and
    purity = N_yes / N  (= the manual validation rate = a direct purity readout).
    This shows: raising D -> fewer cells but higher purity.

(B) N-cost of the gate, using the Ambiguous_BG pool.
    Bergmann in the ML but closer than 45 um were demoted to Ambiguous_BG.
    Lowering the gate would re-admit them. We count how many sit in each distance
    band below 45 um -> how many extra (un-validated, lower-confidence) cells a
    looser gate would add.

Inputs : fc_decisions_zhuang.csv, ../layer_segmentation/BG_assignments_kept.csv
Outputs: threshold_scan.csv, threshold_scan.png
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
dec  = pd.read_csv(os.path.join(HERE, "fc_decisions_zhuang.csv"))
kept = pd.read_csv(os.path.join(HERE, "..", "layer_segmentation", "BG_assignments_kept.csv"))

dec["yes"] = (dec["decision"] == "yes")
N_TOTAL   = len(dec)
N_YES_ALL = int(dec["yes"].sum())
print(f"swipe candidates: {N_TOTAL}  |  human-confirmed yes: {N_YES_ALL}  "
      f"(baseline purity at 45 um = {N_YES_ALL/N_TOTAL:.3f})")

# ---- (A) purity vs distance threshold (gold-label validation rate) ----
D_grid = [45, 50, 55, 60, 65, 70, 80, 90, 100, 120, 150]
rowsA = []
for D in D_grid:
    m = dec["dist_to_PCL_um"] >= D
    n = int(m.sum()); nyes = int((m & dec["yes"]).sum())
    rowsA.append(dict(threshold_um=D, n_candidates=n, n_confirmed=nyes,
                      purity=round(nyes / n, 3) if n else np.nan,
                      pct_of_gold_kept=round(100 * nyes / N_YES_ALL, 1)))
A = pd.DataFrame(rowsA)
print("\n(A) purity vs distance (gold labels):")
print(A.to_string(index=False))

# ---- (B) N-cost of relaxing the 45 um gate (Ambiguous_BG re-admitted) ----
amb = kept[kept["region"] == "Ambiguous_BG"]
bands = [(30, 45), (20, 30), (15, 20), (8, 15)]
rowsB = []
cum = 0
for lo, hi in bands:
    n = int(((amb["dist_to_PCL_um"] >= lo) & (amb["dist_to_PCL_um"] < hi)).sum())
    cum += n
    rowsB.append(dict(band_um=f"[{lo},{hi})", n_ambiguous_readmitted=n,
                      cumulative_extra_if_gate=lo, cumulative_extra_cells=cum))
B = pd.DataFrame(rowsB)
print("\n(B) extra (un-validated) cells if the 45 um gate were lowered:")
print(B.to_string(index=False))

# ---- save combined CSV ----
out_csv = os.path.join(HERE, "threshold_scan.csv")
with open(out_csv, "w") as f:
    f.write("# (A) purity vs distance threshold (gold-label validation rate)\n")
    A.to_csv(f, index=False)
    f.write("\n# (B) extra un-validated cells re-admitted if the 45um gate is lowered\n")
    B.to_csv(f, index=False)
print("\nsaved", out_csv)

# ---- figure ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

# left: purity & N vs threshold
axb = ax1.twinx()
axb.bar(A["threshold_um"], A["n_candidates"], width=3.2, color="#cfd8e3",
        label="N candidats (≥ seuil)", zorder=1)
axb.bar(A["threshold_um"], A["n_confirmed"], width=3.2, color="#9fb4cc",
        label="N confirmés (gold)", zorder=2)
ax1.plot(A["threshold_um"], A["purity"], "-o", color="#d62728", lw=2,
         label="pureté = taux de validation", zorder=3)
ax1.axvline(45, color="k", ls="--", lw=0.8)
ax1.annotate("seuil actuel\n45 µm", (45, 0.5), fontsize=8, ha="left",
             xytext=(50, 0.55), color="k")
ax1.set_xlabel("seuil de distance à la PCL (µm)")
ax1.set_ylabel("pureté (fraction validée à la main)", color="#d62728")
ax1.set_ylim(0.6, 1.02); ax1.tick_params(axis="y", colors="#d62728")
axb.set_ylabel("nombre de cellules", color="#5b6b80")
ax1.set_title("(A) Compromis N ↔ pureté (labels gold)")
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, fontsize=7, loc="lower center")
ax1.set_zorder(axb.get_zorder() + 1); ax1.patch.set_visible(False)

# right: distance distribution, yes vs no
ax2.hist([dec.loc[dec.yes, "dist_to_PCL_um"], dec.loc[~dec.yes, "dist_to_PCL_um"]],
         bins=np.arange(45, 200, 10), stacked=True,
         color=["#2ca02c", "#d62728"], label=["validé (oui)", "rejeté (non)"])
ax2.axvline(45, color="k", ls="--", lw=0.8)
ax2.set_xlabel("dist_to_PCL_um du candidat"); ax2.set_ylabel("candidats")
ax2.set_title("(B) Où tombent les rejets ?"); ax2.legend(fontsize=8)

fig.tight_layout()
out_png = os.path.join(HERE, "threshold_scan.png")
fig.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close()
print("saved", out_png)
