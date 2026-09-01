#!/usr/bin/env python
"""
Article figure palette — SINGLE SOURCE OF TRUTH.
Import this in every figure script so colours are identical everywhere.
Categorical hues validated CVD-safe (dataviz validator, worst adjacent ΔE = 24.2,
white surface). Three fills (Bergmann_oth, Astro_NT, Granule) are <3:1 on white ->
RELIEF RULE: always give them a 0.5 pt dark outline + a label/legend.

    from charte_graphique.palette import GROUP, LAYER, SEMANTIC, seq_blue, div_bwr
"""
from matplotlib.colors import LinearSegmentedColormap, to_rgba

# ---------------------------------------------------------------- entity colours
# ONE colour per group/cell type, used in EVERY figure (schematic, scatter, bars).
GROUP = {
    "FC":           "#d62728",  # red  — Fañanas cell / displaced Bergmann (the focus/accent)
    "BG":           "#2a78d6",  # blue — PCL Bergmann
    "Bergmann_oth": "#1baf7a",  # aqua — reference Bergmann (relief: outline+label)
    "Astro_CB":     "#009E73",  # bluish-green (Okabe-Ito) — cerebellar astrocyte; changed 2026-08-10:
                                # #008300 was near-isoluminant with FC red under deuteranopia (dE~7);
                                # #009E73 keeps dE~41 vs red for deutan viewers (Fig2 audit)
    "Astro_NT":     "#eda100",  # yellow — non-telencephalic astrocyte (relief)
    "Purkinje":     "#4a3aa7",  # violet
    "MLI":          "#eb6834",  # orange — molecular-layer interneuron
    "Granule":      "#e87ba4",  # magenta (relief)
    "Oligo":        "#898781",  # neutral grey — "Other" (9th type, never a hue)
}
# fills that fail the 3:1 contrast floor on white -> must carry a dark outline + label
RELIEF = {"Bergmann_oth", "Astro_NT", "Granule"}
GROUP_EDGE = "#333333"          # standard mark outline

# --------------------------------------------------------------- anatomy / layers
# Muted BACKGROUND tints for layer bands in schematics/overlays (never data series).
LAYER = {
    "ML":  "#e9f1fb",   # molecular layer      (light blue)
    "PCL": "#fdefd8",   # Purkinje cell layer  (cream)
    "GCL": "#efe9f7",   # granule cell layer   (lavender)
    "WM":  "#eceef1",   # white matter         (grey)
    "pia":    "#7f8791",
    "border": "#cfd4da",
}

# ------------------------------------------------------------------- semantic/DE
SEMANTIC = {
    "up_FC":    GROUP["FC"],     # significant, higher in FC
    "up_BG":    GROUP["BG"],     # significant, higher in BG
    "contam":   "#e8a13a",       # contamination flag — ALWAYS with an 'x' marker (2nd encoding)
    "ns":       "#b8b6b0",       # non-significant / background points
    "equiv":    "#5a5a57",       # equivalence band / null reference
}

# --------------------------------------------------------------------- ink/chrome
INK = {
    "primary":   "#0b0b0b",
    "secondary": "#52514e",
    "muted":     "#898781",   # axis labels, ticks
    "grid":      "#e1e0d9",   # hairline gridline
    "axis":      "#c3c2b7",   # baseline / spine
    "surface":   "#ffffff",   # print figure surface (white)
}

# ------------------------------------------------------ sequential ramp (blue)
_SEQ_BLUE = ["#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7",
             "#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281","#0d366b"]
seq_blue = LinearSegmentedColormap.from_list("seq_blue", _SEQ_BLUE)
# second sequential context (aqua), if two heatmaps co-occur
_SEQ_AQUA = ["#d4f2e6","#a9e3cd","#71cfac","#3fba8c","#1baf7a","#159068","#0e6f50"]
seq_aqua = LinearSegmentedColormap.from_list("seq_aqua", _SEQ_AQUA)

# ------------------------------------------------------ diverging (blue↔grey↔red)
# for log2FC etc.: BG(blue) negative pole, neutral grey mid, FC(red) positive pole
div_bwr = LinearSegmentedColormap.from_list(
    "div_bwr", [GROUP["BG"], "#8fb4e6", "#f0efec", "#e79596", GROUP["FC"]])

def rgba(name_or_hex, alpha=1.0):
    return to_rgba(GROUP.get(name_or_hex, name_or_hex), alpha)

# fixed legend/plot ORDER (CVD-optimal; keeps similar hues non-adjacent)
ORDER = ["BG","Bergmann_oth","Astro_NT","Astro_CB","Purkinje","FC","Granule","MLI","Oligo"]

if __name__ == "__main__":
    for k,v in GROUP.items(): print(f"{k:14s} {v}  {'(relief)' if k in RELIEF else ''}")
