# Hao 2024 — Stereo-seq

Dataset spatial principal. Il n’existe qu’un seul état actif dans ce dossier.

## Organisation

| Dossier | Rôle |
|---|---|
| `data_raw/` | 59 sections Stereo-seq H5AD historiques et atlas Kozareva normalisé |
| `analysis/` | Pipeline global, contrôles et objets traités actifs |
| `analysis/depth_matched/` | Analyse FC–BG de référence pour l’article |
| `layer_segmentation/` | Attribution anatomique FC/BG |
| `swipe_validation/` | Décisions manuelles sur les candidats |

L’analyse de référence est lancée avec :

```bash
python3 datasets/hao_2024_stereo_seq/analysis/depth_matched/run_depth_matched_de.py
```

Elle remplace directement `DE_FC_vs_BG.csv`, `module_scores.csv`, `cell_sets.json`,
`zhuang_comparison.csv` et `volcano.png`. Aucun fichier numéroté ou copie de version ne doit
être créé.

## Résultat courant

- La cohorte active est **Mouse1 uniquement** : 89 FC confirmées et **7 230 BG confiantes** dans le
  recensement anatomique complet. La cohorte d’expression contient 89 FC, un sous-échantillon de
  1 524 BG et 960 astrocytes parenchymateux ; elle ne doit pas servir à estimer la fréquence anatomique.
- Les distances physiques sont calibrées sur la grille 10 µm de Mouse1 (médiane FC 86 µm, n=89).
  Mouse2 emploie une autre échelle de coordonnées non calibrée : elle est exclue de toute cohorte,
  table, statistique et figure active (critère géométrique, aveugle à l'expression — la définition FC
  est métrique ; justification complète au canonique, review point 34). Ses H5AD et décisions
  manuelles restent conservés comme provenance, ainsi que le contrôle directionnel
  `scripts/05_mouse2_validation.py` / `results/cross_mouse_effect_sizes.csv` : concordant avec
  Mouse1 (85 % de signes concordants à |d_M1|≥0,2, r=0,52) — l'exclusion n'est pas protectrice
  du résultat.
- Analyse primaire d’expression : Mouse1, profondeur égalisée à 600 transcrits par cellule,
  permutation intra-section.
- Identité Bergmann conservée, aucun type FC distinct et lean astrocytaire nul à profondeur égale.

Pour tout chiffre inter-datasets, utiliser uniquement
`../../comparison/CANONICAL_NUMBERS.md`.

## Article source

`Hao_2024_Science_cerebellar_atlas.pdf` — Hao S. *et al.* (2024) *Cross-species single-cell spatial
transcriptomic atlases of the cerebellar cortex.* **Science** 385:eado3927, doi 10.1126/science.ado3927.
Déposé le 2026-08-21 via l'accès institutionnel (l'article n'est pas en accès libre : ni PMC, ni
préprint bioRxiv). 17 pages — Research Article Summary + texte intégral.
