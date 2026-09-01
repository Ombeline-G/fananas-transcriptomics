# Suyama 2026 — Xenium adulte FC versus BG

Dataset spatial principal Xenium Prime 5K. L'article source traite P0/P6 ; le projet utilise
la section adulte supplémentaire déposée avec l'archive Zenodo.

## Données conservées

- `data/P0_P6_Adult_mousecerebellum_5K.tar.gz` : archive brute Zenodo, vérifiée par `gzip -t`
  le 12 août 2026.
- `data/full_slide.rds` et `data/adult_cerebellum.rds` : objets nécessaires à la reconstruction
  de la section adulte.
- `data/adult_laminar.rds` : géométrie laminaire adulte.
- `results/tables/` : exports nécessaires à l'harmonisation.
- `results/swipe_validation/` et `results/fc_decisions_suyama.csv` : validation manuelle.

La copie brute extraite de 28 Go a été supprimée parce qu'elle dupliquait l'archive vérifiée.
Le pipeline local courant repart de `full_slide.rds`; l'archive brute reste la source froide à
extraire uniquement si cet objet doit être reconstruit.

- `data/cell_feature_matrix/` (features/barcodes/matrix.mtx.gz, ~600 Mo) : ré-extrait de
  l'archive Zenodo le 17 août 2026 pour `comparison/figure_scripts/build_equivalence_scoping.py`
  (review point 7), qui recalcule log2(1+CP10k) sur les 5 006 gènes des 3 672 astroglies
  canoniques sans passer par R/Seurat (validé exactement contre `adult_panel_expr.csv` sur les
  61 gènes partagés ; totaux par cellule = features « Gene Expression » seules).

Les objets P0/P6 dérivés et les scripts développementaux 01–05 ont été retirés : ils ne
participent pas à la comparaison adulte.

## Pipeline adulte

- `scripts/08` à `11` : crop adulte, clustering, géométrie laminaire et candidats initiaux.
- `scripts/17`, `19`, `21` à `23` : sélection des candidats et support de swipe.
- `scripts/30` à `32` : exports d'expression nécessaires aux analyses harmonisées.
- `comparison/figure_scripts/build_suyama_purity_gate.py` : construction du set FC canonique
  (pureté + anti-PCL ≥ 30 µm) et du contrôle strict 45 µm.
- `comparison/figure_scripts/suyama_build.py` : parquets harmonisés finaux.

Les scripts R détectent la racine lorsqu'ils sont lancés avec `Rscript`. En mode interactif,
définir `BGROOT` vers la racine du projet.

## Résultat autorisé

- 73 candidates confirmées au swipe, 56 après gate de pureté symétrique pré-spécifié, puis **44 FC
  canoniques** après la règle anti-PCL (soma à d_pkj ≥ 30 µm du soma de Purkinje le plus proche).
- 1 823 BG et 1 805 astrocytes parenchymateux.
- Identité Bergmann conservée ; aucun type FC distinct.
- Lean astrocytaire harmonisé +0,24 sur 42 paires, p=8,4e-4.
- Contrôle strict 45 µm : parmi les 44 FC canoniques, **20** sont à ≥45 µm
  (`results/tables/suyama_fc_gate_strict45.json`). Sur ces 20 paires appariées, le lean
  garde sa direction mais perd la significativité (Δ module +0,12, p=0,39 ; Aqp4 +0,24,
  p=0,08 ; `../../comparison/tables/suyama_strict45_summary.csv`).

La source de chiffres est `../../comparison/CANONICAL_NUMBERS.md`.

## Article source

`Suyama_2026_Glia_bergmann_glia_like_progenitors.pdf` — Suyama K. *et al.* (2026) *Molecular
Characteristics and Differentiation Control Mechanisms of Bergmann Glia-Like Progenitors in the
Postnatal Mouse Cerebellum.* **Glia** 74:e70160, doi 10.1002/glia.70160. 18 pages.
Données : dépôt Zenodo 16562843 (CC-BY 4.0). Renommé le 2026-08-21 (le nom précédent était un
reliquat de téléchargement tronqué ; vérifié : aucun script ne le référençait).
