# Validation manuelle FC — outil « swipe » (source : segmentation)

Curation manuelle (gold standard) des **FC** définis par la segmentation géométrique
(`../layer_segmentation/`), sur les 20 coupes Zhang validées.

## Candidats
- **240 FC** = Bergmann-subclass dont le soma est clairement dans la **ML**,
  **≥45 µm de la PCL** (clairement déplacé), après **déclumping** (les amas retirés).
- Triés par **confiance** (segmentation) décroissante.

## Ouvrir
- Double-clic sur `swipe.html` (données embarquées dans `swipe_zhuang.js`, hors-ligne), ou via le preview Claude Code (port 8765).

## Lire la carte
- **★ rouge** = soma candidat FC (centré). Fond = **masques de segmentation** : ML (bleu), PCL (orange), GCL (violet), WM (gris). Points bleu foncé = BG ; orange = Purkinje.
- Coupe entière (gauche) + zoom local (droite). Panneau droit : confiance, distance à la PCL, support Purkinje.

## Décider
✓ Oui (→/Y) · ? Peut-être (↑/M) · ✗ Non (←/N) · Retour = annuler. Sauvegarde auto.
**Critère :** ★ bien dans la ML (bleu), au-dessus de la ligne de Purkinje (orange), isolé (pas dans un amas).

## Export
Bouton **« Exporter les décisions (CSV) »** → `fc_decisions_zhuang.csv`
(colonnes : cell_label, section, animal, confidence, dist_to_PCL_um, decision).
Renvoie-le-moi → set FC « gold » → convergence Hao↔Zhuang sur cellules validées.

## Régénérer (si besoin)
`python gen_swipe_from_segmentation.py` (lit `../layer_segmentation/BG_assignments_kept.csv` + masks).

## ABCA-1 (review point 35)

Quatrième animal, option A (valider et inclure). Ouvrir **`swipe_abca1.html`** (double-clic,
hors-ligne) : **406 candidats** sur 14 sections, mêmes critères, même interface, même aveuglement
(géométrie seule). Décisions stockées séparément (`bgml_swipe_decisions_v3seg_abca1`) ; export →
`bgml_decisions_abca1.csv`. Données régénérables : `python gen_swipe_from_segmentation.py abca1`
(lit `../layer_segmentation/FC_candidates_ABCA1_for_swipe.csv` + `BG_assignments_ABCA1.csv`).
Après swipe : intégration ABCA-1 par le pipeline identique (voir CANONICAL, extension
« Intégration ABCA-1 post-swipe »).

## Bande 35–45 µm (review S. Alon, novembre 2026)

Analyse de sensibilité à la marge de sélection : les Bergmann du masque ML situées entre 35 et 45 µm
de la PCL n'avaient jamais été relues (la marge de l'étude est 45 µm). Ouvrir **`swipe_band35.html`**
(double-clic, hors-ligne) : **352 candidats** sur 33 sections des quatre animaux (ABCA-1 : 236 ·
ABCA-2 : 4 · ABCA-3 : 94 · ABCA-4 : 18), mêmes critères, même interface, même aveuglement (géométrie
seule), tri par confiance décroissante. Décisions stockées séparément
(`bgml_swipe_decisions_v3seg_band35`) ; export → **`fc_decisions_band35.csv`**, à déposer dans ce
dossier, toutes décisions résolues en oui/non. Données régénérables :
`python gen_swipe_from_segmentation.py band35` (écrit aussi
`../layer_segmentation/FC_candidates_band35_45_for_swipe.csv`).

**Swipe fait le 2026-11-23 (O.G.) : 74 acceptées / 352** (ABCA-1 26/236 · ABCA-2 1/4 · ABCA-3 37/94 ·
ABCA-4 10/18), décisions dans `fc_decisions_band35.csv`. Ces 74 cellules n'entrent que dans la cohorte
de sensibilité à 35 µm (`comparison/figure_scripts/build_selection_sensitivity.py`, Figure S10) ; la
cohorte de l'étude reste celle des 362 FC à ≥45 µm (canonique, addendum 60).
