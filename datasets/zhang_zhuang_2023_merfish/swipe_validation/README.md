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
