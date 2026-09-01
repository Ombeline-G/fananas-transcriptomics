# Zhang/Zhuang 2023 — MERFISH FC versus BG

Dataset spatial principal issu de l'Allen Brain Cell Atlas. La citation correcte est
**Zhang et al. 2023, Nature** ; « Zhuang » désigne le laboratoire et le nom du dataset.

## Pipeline actif

1. `07_pool_layers_all4.py` — métadonnées cérébelleuses mutualisées.
2. `08_download_expr_2to4.py` — matrices d'expression MERFISH.
3. `14_multivote.py` — préparation des quatre critères anatomiques.
4. `layer_segmentation/` — segmentation GCL-anchored, seuil 45 µm et QC.
5. `swipe_validation/` — 190 FC retenues sur 240 candidates.
6. `validated_analysis/` — DE validée.
7. `identity_analysis/` — tests positifs de conservation Bergmann.
8. `hierarchical_stats/` — inférence par coupe/animal et correction de pseudoréplication.
9. `ml_glia_inventory/` — census sans filtre Bergmann préalable.

## Résultat autorisé

- 190 FC, 3 570 BG, 1 427 astrocytes parenchymateux ; 3 animaux, 20 sections.
- Identité Bergmann conservée ; absence de type FC distinct.
- Léger lean astrocytaire harmonisé : +0,19 sur 182 paires, p=1,4e-5.
- ABCA-1 : 406 candidates pré-swipe ; résultat exploratoire non citable avant validation.
  **Outil de validation prêt** (review point 35) : `swipe_validation/swipe_abca1.html` (double-clic,
  hors-ligne ; 406 cartes, 14 sections, export `bgml_decisions_abca1.csv`) ; régénérable par
  `python gen_swipe_from_segmentation.py abca1`. Après swipe : intégration par le pipeline identique.

La source de chiffres inter-datasets est `../../comparison/CANONICAL_NUMBERS.md`.

## Article source

`Zhang_2023_Nature_whole_mouse_brain_atlas.pdf` — Zhang M. *et al.* (2023) *Molecularly defined and
spatially resolved cell atlas of the whole mouse brain.* **Nature** 624:343–354, doi 10.1038/s41586-023-06808-9.
Téléchargé le 2026-08-21 depuis nature.com (accès libre, **CC BY**). Rappel de citation : **Zhang et al.**,
jamais « Zhuang et al. » (Zhuang = laboratoire).
