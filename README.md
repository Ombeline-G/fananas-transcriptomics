# Fañanas cells of the adult mouse cerebellum — analysis code

Analysis code, computed tables, and manual-validation materials for the transcriptomic
comparison of Fañanas cells (FC) and Bergmann glia (BG) across three spatial transcriptomic
datasets (MERFISH, Stereo-seq, Xenium) and two dissociated atlases (Kozareva 2021, Siletti 2023).

The associated manuscript reports that FC are a rare Bergmann-class astroglial population with
somata displaced into the molecular layer, not a distinct transcriptomic type; a modest,
position-associated astrocytic expression shift is reported as a secondary finding.

## Layout

| Path | Content |
|---|---|
| `comparison/figure_scripts/` | All analysis and figure scripts (Python; entry order in `comparison/README.md`) |
| `comparison/tables/` | Computed tables — outputs and the stable inputs the rebuild leaves in place |
| `comparison/rescore_package/` | Second-rater validation: mini swipe tools, both raters' decision files, official merge, Cohen's κ pipeline |
| `comparison/requirements.txt` | Pinned production environment (Python 3.13.7) |
| `charte_graphique/` | Figure style single source of truth (palette, mplstyle), imported by `figstyle.py` |
| `datasets/*/swipe_validation/` | Interactive swipe tools, generators, and per-candidate manual decisions (1,152 candidates) |
| `datasets/*/README.md` | Per-dataset local pipeline documentation and data provenance |

## Reproducing

```bash
pip install -r comparison/requirements.txt
export BGROOT="/path/to/this/repository"   # optional; paths also resolve relative to the scripts
```

Then run the table builders and figure scripts in the order given in `comparison/README.md`
(dependency constraints are documented there).

Two deliberate environment notes, declared for reproducibility:

- **`umap-learn` is not part of the production environment.** `Figure3_identity.py` then falls
  back to plotting the first two principal components; the published panel 3d **is** that PCA
  projection, as its axis labels state.
- **`scanpy` is not part of the production environment.** `build_atlas_siletti.py` is therefore
  not re-executable as-is; its outputs (`atlas_*_verdict.csv`, …) ship in `comparison/tables/`
  as stable inputs.

## Manual validation (swipe)

Every geometric FC candidate was accepted or rejected in a dedicated offline interface
(`datasets/*/swipe_validation/swipe*.html`; geometry only, blinded to expression). All
per-candidate decisions are included (`fc_decisions_*.csv`). `comparison/rescore_package/`
contains the second-rater re-scoring of a seed-fixed, stratified 15% sample (173 candidates):
the four restricted tools, both raters' decisions, the official merge (`rater2_all.csv`), and
`comparison/figure_scripts/compute_rescore_kappa.py`, which reproduces the inter-rater
agreement reported in Methods (pooled Cohen's κ = 0.46).

## Data

Raw datasets are not redistributed here; each `datasets/*/README.md` documents its source.
The Xenium dataset (Suyama et al. 2026) is on Zenodo (record 16562843, CC-BY 4.0); MERFISH is
from the Allen Brain Cell Atlas (Zhang et al. 2023); Stereo-seq from Hao et al. 2024; the
dissociated atlases are Kozareva et al. 2021 (GSE165371) and Siletti et al. 2023.

The manuscript, internal decision logs, and figure outputs are not part of this repository
(references to them in `comparison/README.md` point to the project working tree).

## License

Code is released under the MIT License (see `LICENSE`). The datasets remain under their
original licenses.
