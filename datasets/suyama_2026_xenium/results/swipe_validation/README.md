# Validation manuelle FC — outil « swipe » (Suyama adulte)

Curation manuelle (gold standard) des **145 FC** définis par la segmentation géométrique
(`d_gran` = distance à la couche granulaire ≥ 30 µm), sur le **cervelet adulte** Suyama (1 coupe).

## Candidats
- **145 FC** = Bergmann dont le soma est dans la **ML** (≥30 µm de la couche granulaire).
- Triés par **confiance** (`d_gran`) décroissante : carte 1 = plus franchement ML (~96 µm),
  carte 145 = limite (~30 µm).

## Ouvrir
- Double-clic sur `swipe.html` (données embarquées, hors-ligne), ou via le preview Claude Code (port 8767).

## Lire la carte
- **★ rouge** = soma candidat FC (centré). Fond = couches : ML (bleu), GCL (violet), WM (gris).
  Points orange = Purkinje (= la PCL) ; bleu foncé = BG_PCL.
- Cervelet entier (gauche, rectangle = zone du zoom) + zoom local 200 µm (droite).
- Panneau droit : distance à la granule, distance au Purkinje, nb de Purkinje < 40 µm.

## Décider
✓ Oui (→/Y) · ? Peut-être (↑/M) · ✗ Non (←/N) · Retour = annuler. Sauvegarde auto.
**Critère :** ★ bien dans la ML (bleu), au-dessus de la ligne de Purkinje (orange), isolé.

## Export
Bouton **« Exporter les décisions (CSV) »** → `fc_decisions_suyama.csv`
(colonnes : cell_label, dataset, section, animal, d_gran_um, d_purkinje_um, decision).
Renvoie-le-moi → set FC « gold » → je relance contrôle soma + DE + tendance sur les cellules validées.

## Régénérer (si besoin)
`Rscript ../../scripts/22_swipe_data.R` puis `python ../../scripts/23_build_swipe.py`.
