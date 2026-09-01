# Comparaison active

Ce dossier contient l’unique analyse inter-datasets du projet.

## Sources de vérité

1. `CANONICAL_NUMBERS.md` — chiffres et conclusions.
2. `tables/` — données calculées soutenant ces chiffres.
3. `notes/HARMONIZATION_protocol.md` — protocole commun.
4. `notes/DECISIONS.md` — décisions éditoriales actives.
5. `notes/FIGURE_PLAN.md` et `notes/figure_qc/` — figures et validations.

## Sorties finales

`figures/` contient les figures actives : Figure1–5 et les figures supplémentaires S1–S6 et S8–S16
(S7 fusionnée dans S6 au point 27 de la review ; renumérotation à l'assemblage). Chaque figure existe
en PDF, SVG et PNG. Aucun autre dossier de figures ne fait foi.

| Figure principale | Message | Statut |
|---|---|---|
| Figure 1 | Localisation des somas FC | Validée après audit des sources anatomiques |
| Figure 2 | Census position-first de la couche moléculaire | Validée |
| Figure 3 | Identité Bergmann conservée | Validée |
| Figure 4 | Absence de type FC distinct | Validée |
| Figure 5 | Identité conservée, état sensible à la position (lean astrocytaire) | Validée |

## Reconstruction

Depuis la racine du projet, d’abord les builders de tables, puis les figures :

```bash
python3 comparison/figure_scripts/build_suyama_purity_gate.py
python3 comparison/figure_scripts/suyama_build.py
python3 comparison/figure_scripts/hao_build.py
python3 comparison/figure_scripts/hao_func_export.py
python3 comparison/figure_scripts/restrict_hao_census_mouse1.py
python3 comparison/figure_scripts/harmonized_lean.py
python3 comparison/figure_scripts/build_identity_metrics.py
python3 comparison/figure_scripts/knn_sensitivity.py
python3 comparison/figure_scripts/build_full_de.py
python3 comparison/figure_scripts/build_de_canonical.py
python3 comparison/figure_scripts/build_spillover.py
python3 comparison/figure_scripts/build_equivalence.py
python3 comparison/figure_scripts/build_equivalence_scoping.py
python3 comparison/figure_scripts/build_normalized_depth.py
python3 comparison/figure_scripts/build_fc_by_lobule.py
python3 comparison/figure_scripts/build_lean_depth_continuum.py
python3 comparison/figure_scripts/build_vascular_proximity.py
python3 comparison/figure_scripts/build_annotation_free_identity.py
python3 comparison/figure_scripts/build_dotplot.py
python3 comparison/figure_scripts/build_module_summary.py
python3 comparison/figure_scripts/build_atlas_kozareva.py
python3 comparison/figure_scripts/build_atlas_siletti.py
python3 comparison/figure_scripts/Figure1_localization.py
python3 comparison/figure_scripts/Figure2_ML_census.py
python3 comparison/figure_scripts/FigureS4_signatures.py
python3 comparison/figure_scripts/Figure3_identity.py
python3 comparison/figure_scripts/FigureS3_confounds.py
python3 comparison/figure_scripts/Figure4_nottype.py
python3 comparison/figure_scripts/Figure5_lean.py
python3 comparison/figure_scripts/FigureS1_layers.py
python3 comparison/figure_scripts/FigureS2_census_hao.py
python3 comparison/figure_scripts/FigureS5_heterogeneity.py
python3 comparison/figure_scripts/FigureS6_S8_S10.py
python3 comparison/figure_scripts/FigureS9_regression.py
python3 comparison/figure_scripts/FigureS11_annotation_free.py
python3 comparison/figure_scripts/FigureS12_atlases.py
python3 comparison/figure_scripts/FigureS13_lobules.py
python3 comparison/figure_scripts/FigureS14_continuum.py
python3 comparison/figure_scripts/FigureS15_vessels.py
python3 comparison/figure_scripts/FigureS16_swipe_examples.py
```

Contraintes d’ordre vérifiées dans les scripts : `FigureS4_signatures.py` s’exécute avant
`Figure3_identity.py` (qui lit `fig6_all_contrasts_DE.csv`, produit par S4) et avant
`FigureS3_confounds.py` (qui ajoute une colonne à `Figure6_signature_table.csv`, produit par S4) ;
chaque builder précède la figure qui lit ses tables (S6/S8/S10 ← `build_module_summary.py`,
S11 ← `build_annotation_free_identity.py`, S12 ← `build_atlas_kozareva.py` +
`build_atlas_siletti.py`, S13 ← `build_fc_by_lobule.py`, S14 ← `build_lean_depth_continuum.py`,
S15 ← `build_vascular_proximity.py`) ; `build_full_de.py` consomme le gate et les paires
appariées écrits par `build_suyama_purity_gate.py` ; `build_dotplot.py` consomme
`hao_bg_gene_stats.csv` écrit par `hao_func_export.py`.

`zhuang_build.py` n’est nécessaire que si la cohorte Zhuang change. Hao est explicitement reconstruit
ici afin de garantir que la cohorte active reste limitée à Mouse1. Certaines tables de `tables/` sont
des entrées stables sans builder dans ce dossier et la reconstruction les laisse en place : census
position-first Zhuang (`zhuang_ML_census_cells.parquet`, `zhuang_census_identity*.csv`), census Hao
(`hao_ML_census_expr.parquet`, restreint en place par `restrict_hao_census_mouse1.py`),
`robust_layers_zhuang.parquet`, `hao_fc_share_by_section.csv`, `fig7_functional_scores.csv`,
`threshold_robustness.csv` et les verdicts atlas (`atlas_siletti_donor_verdict.csv`,
`atlas_kozareva_sex_verdict.csv`).

`build_swipe_rescore.py` et `compute_rescore_kappa.py` sont l’infrastructure de double lecture de la
validation swipe (κ de Cohen, review point 26) : le premier régénère l’échantillon à re-scorer
(`swipe_rescore_sample.csv`), le second attend en argument le fichier du second lecteur ; ils ne font
pas partie de la reconstruction standard.

Le manuscrit courant est `MANUSCRIPT.md`, construit depuis `CANONICAL_NUMBERS.md` et
`notes/FIGURE_CAPTIONS_draft.md` ; les reviews sont consignées dans `REVIEW_LOG_Antonin_v1.md` et
`REVIEW_LOG_Ombeline.md`. Le manuscrit ne contient plus de section interne « Open items » : le
suivi (dont le bloqueur de soumission κ second lecteur) vit dans `notes/OPEN_ITEMS.md`.

La copie de lecture `MANUSCRIPT_print.pdf` se régénère depuis `MANUSCRIPT.md` avec :

```bash
/Users/ombelinegruere/.python-global/bin/python comparison/figure_scripts/md2pdf.py
```

(reportlab + Pillow ; polices Arial du système. Ce script rend un sous-ensemble de Markdown —
titres, gras/italique, listes, tables, images — et intègre les figures depuis `figures/`.)
