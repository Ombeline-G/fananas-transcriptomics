#!/usr/bin/env python
"""Suyama non-cortical mask — shared by all Figure-1 panels.

The adult Xenium crop (Suyama 2026) includes a left + bottom tissue region that is NOT cerebellar
cortex: it contains no granule mass and no Purkinje cells (the two hallmark cortical elements),
being 'Other neurons' + astrocytes/oligo/vascular — adjacent non-cerebellar tissue (brainstem /
deep nuclei) caught in the crop. It must not be interpreted, so Figure-1 panels exclude it.

Definition (purely spatial, reproducible): build the granule-cell MASS as a 2-D region (bin the
section, local granule fraction, smooth, threshold, morphology) and drop every cell whose distance
to that mass exceeds CORTEX_MARGIN_UM. Distance to the granule *mass* (not to the nearest granule
*cell*) is the right metric — a few stray/misclassified granule near the block would otherwise let
the block escape a nearest-cell threshold. The molecular layer reaches ~175 µm from the mass, so a
200 µm margin keeps the whole cortex ribbon (GCL + PCL + ML + adjacent WM) and cuts only beyond it.

Verified 2026-07-15: removes 6 200 cells (Other 3 600, Astrocyte 1 187, Oligo 672, Vascular 457,
OPC 146, Immune 116, + 19 'granule'/3 'MLI' that are MIS-ANNOTATED — 0 granule neighbours within
40 µm, i.e. not real cortical cells). PROTECT keeps Bergmann + Purkinje so **all Fañanas cells
(Bergmann in the ML) are retained** (0 Bergmann / 0 Purkinje ever dropped). A small intra-folia
vascular cluster (a legitimate cerebellar vessel, dm small) is kept. Suyama typed 40 547 -> 34 347.
"""
import json
import numpy as np
import pandas as pd
from scipy.ndimage import (gaussian_filter, binary_closing, binary_opening,
                           binary_fill_holes, label as ndlabel, distance_transform_edt)

BIN = 30.0                 # grid resolution (µm)
SMOOTH = 40.0              # gaussian smoothing of the granule-fraction field (µm)
GRAN_THR = 0.30            # local granule fraction to call a bin 'granule mass'
MIN_REGION_UM2 = 13000.0   # drop granule-mass islands smaller than this
CORTEX_MARGIN_UM = 200.0   # cells farther than this from the granule mass = non-cortical
PROTECT = ("Bergmann", "Purkinje")   # PCL soma types never dropped -> guarantees the FC are kept


def _granule_mass_distance(d):
    """Per-cell distance (µm) to the granule-cell mass region."""
    X = d.x.values; Y = d.y.values
    ix = ((X - X.min()) / BIN).astype(int); iy = ((Y - Y.min()) / BIN).astype(int)
    nx, ny = ix.max() + 1, iy.max() + 1
    def grid(v):
        G = np.zeros((nx, ny)); np.add.at(G, (ix, iy), v); return G
    sig = SMOOTH / BIN
    tot = gaussian_filter(grid(np.ones(len(X))), sig)
    fracG = gaussian_filter(grid((d.ct.values == "Granule").astype(float)), sig) / (tot + 1e-6)
    gm = fracG > GRAN_THR
    gm = binary_closing(gm, iterations=2); gm = binary_opening(gm, iterations=2); gm = binary_fill_holes(gm)
    lab, _ = ndlabel(gm); sz = np.bincount(lab.ravel()); sz[0] = 0
    gm = np.isin(lab, np.where(sz >= int(MIN_REGION_UM2 / (BIN * BIN)))[0])
    return (distance_transform_edt(~gm) * BIN)[ix, iy]


def noncortical_cell_ids(base):
    """Return the set of Suyama cell IDs (as str) that fall outside the cerebellar cortex."""
    su = f"{base}/datasets/suyama_2026_xenium"
    d = pd.read_csv(f"{su}/results/tables/adult_cells_scores.csv")[["cell", "x", "y", "cl"]]
    d["ct"] = d.cl.map(json.load(open(f"{su}/results/tables/adult_cl2type.json")))
    d = d.dropna(subset=["ct"]).reset_index(drop=True)
    dm = _granule_mass_distance(d)
    noncort = (dm > CORTEX_MARGIN_UM) & (~d.ct.isin(PROTECT).values)
    return set(d.loc[noncort, "cell"].astype(str))
