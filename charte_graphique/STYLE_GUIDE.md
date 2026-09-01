# Charte graphique — figures de l'article « FC = Bergmann déplacées »

Spec **complète et explicite** : polices, tailles (pt), couleurs (hex + RGB), épaisseurs (pt).
**Toute valeur est fixée ici** — aucune n'est laissée à l'interprétation.
**Emplacement :** `8_Transcriptomic_comparison_datasets/charte_graphique/` — à la racine du dossier de comparaison, commune à **tous** les jeux (Zhuang, Hao2024, Suyama2026, comparison).
Source de vérité machine : [`palette.py`](palette.py) + [`article.mplstyle`](article.mplstyle).
Règle absolue : **aucune couleur ni taille écrite en dur** dans les scripts.

```python
# depuis un script dans un sous-dossier de jeu de données (…/spatial_transcriptomic_*/)
import os, sys, matplotlib.pyplot as plt
CHARTE = os.path.join(os.path.dirname(__file__), "..", "charte_graphique")
plt.style.use(os.path.join(CHARTE, "article.mplstyle"))
sys.path.insert(0, os.path.dirname(CHARTE))
from charte_graphique.palette import GROUP, LAYER, SEMANTIC, INK, seq_blue, seq_aqua, div_bwr, ORDER, RELIEF, GROUP_EDGE
```

---

## 0. Principe directeur — LISIBILITÉ À PETITE TAILLE (prioritaire sur tout le reste)
La figure doit rester **entièrement lisible une fois réduite** (colonne étroite, diapo, vignette).
En cas de conflit avec toute autre règle ci-dessous, **ce principe gagne** : on agrandit ou on retire du contenu, **jamais on ne rétrécit**.
1. **Aucun petit détail.** Aucun élément porteur d'information ne doit disparaître à la réduction. Ce qui compte est **gros, épais, contrasté**.
2. **Texte gros.** Rien < **7 pt** à taille finale ; viser **8–9 pt** pour axes/labels, **12 pt** pour lettres de panneaux. Préférer trop gros à trop petit.
3. **Traits épais.** Aucune ligne < **0,75 pt** ; lignes de données **2–2,5 pt** ; marqueurs **≥ 4 pt**.
4. **Peu d'éléments par panneau.** Panneau chargé → le **scinder**, pas le rétrécir. Jamais 6 panneaux entassés dans une figure.
5. **Labels directs > légendes minuscules** ; **fort contraste** ; **espace blanc généreux** (ne pas remplir chaque pixel).
6. **Test de réduction — OBLIGATOIRE avant de valider toute figure :** réduire à **50 %** (ou afficher en ~6 cm de large). Si un texte, un trait ou un marqueur devient difficile à lire → **ÉCHEC** → agrandir/simplifier. Seuils : rien < 7 pt, aucun trait < 0,75 pt, aucun marqueur < 4 pt à taille finale.

---

## 1. Format & export
| Paramètre | Valeur exacte |
|---|---|
| Largeur 1 colonne | **89 mm** = 3,50 in = 1 050 px @300 dpi |
| Largeur 1,5 colonne | **120 mm** = 4,72 in |
| Largeur 2 colonnes | **183 mm** = 7,20 in = 2 160 px @300 dpi |
| Hauteur maximale | **247 mm** = 9,72 in |
| Échelle de dessin | **taille finale 1:1** (jamais de mise à l'échelle après export) |
| Fond figure | blanc `#ffffff` (opaque) |
| Format maître | **PDF** vectoriel, polices **non tracées** (`pdf.fonttype = 42`) |
| Format édition | **SVG** (`svg.fonttype = none` → texte éditable) |
| Format raster | **PNG 600 dpi** (`savefig.dpi = 600`) |
| Rognage | `bbox_inches = tight`, `pad_inches = 0.02` |
| Espace colorimétrique | RGB à l'écran ; **convertir en CMYK** à la soumission si le journal l'exige |
| Gouttière inter-panneaux | **4 mm** |
| Marge de données interne | 2 % (`axes.xmargin = axes.ymargin = 0.02`) |

---

## 2. Polices
| Élément | Police | Graisse | Style |
|---|---|---|---|
| Tout texte | **Arial** (repli : Helvetica → Liberation Sans → DejaVu Sans) | selon §3 | roman |
| Symboles de gènes | Arial | normal | **italique** (*Gfap*, *Agt*, *Calb1*) |
| Variables stat | Arial | normal | **italique** (*P*, *n*, *r*, *t*, *q*) |
| Lettres de panneaux | Arial | **gras (700)** | roman, minuscule |
| Nombres alignés (axes, tableaux) | Arial | normal | `tabular-nums` |

- **Une seule famille** partout. **Deux graisses seulement** : normal (400) et gras (700). Aucune autre.
- Aucune serif, aucune police display, aucune small-caps.

---

## 3. Tailles de texte (pt, à taille finale — **minimum absolu 7 pt**, cf. §0)
Barème « lisibilité d'abord » : texte volontairement **gros** pour survivre à la réduction.
| # | Élément | Taille | Graisse |
|---|---|---|---|
| 1 | Lettre de panneau (a, b, c) | **12** | gras |
| 2 | Titre de figure (si présent) | **12** | gras |
| 3 | Titre de panneau (souvent omis) | **10** | gras |
| 4 | Titre d'axe (x et y) | **9** | normal |
| 5 | Graduations (tick labels) | **8** | normal |
| 6 | Légende | **8** | normal |
| 7 | Étiquette de colorbar | **9** | normal |
| 8 | Graduations de colorbar | **8** | normal |
| 9 | Annotation / label de gène | **8** | normal (gène italique) |
| 10 | Valeur de P / n / étoiles | **8** | normal |
| 11 | Label de barre d'échelle | **8** | normal |
| 12 | Note bas de figure (« not to scale ») | **7** | italique |

---

## 4. Couleurs (hex + RGB)

### 4a. Entités — catégoriel (validé CVD, ΔE 24,2)
| Entité | Hex | RGB | Teinte | Relief¹ |
|---|---|---|---|---|
| **FC** (Fañanas / Bergmann déplacée) | `#d62728` | 214, 39, 40 | rouge | — |
| **BG** (Bergmann PCL) | `#2a78d6` | 42, 120, 214 | bleu | — |
| **Bergmann_oth** (réf.) | `#1baf7a` | 27, 175, 122 | aqua | **oui** |
| **Astro_CB** (contraste clé) | `#008300` | 0, 131, 0 | vert | — |
| **Astro_NT** | `#eda100` | 237, 161, 0 | jaune | **oui** |
| **Purkinje** | `#4a3aa7` | 74, 58, 167 | violet | — |
| **MLI** | `#eb6834` | 235, 104, 52 | orange | — |
| **Granule** | `#e87ba4` | 232, 123, 164 | magenta | **oui** |
| **Oligo** = « Autre » | `#898781` | 137, 135, 129 | gris | — |
| Contour de marque | `#333333` | 51, 51, 51 | — | — |

¹ **Relief** : fill < 3:1 sur blanc → **toujours** contour `#333` 0,5 pt **+** label/légende. Jamais couleur seule.
**Ordre fixe** (légende/tracé, CVD-optimal) : BG · Bergmann_oth · Astro_NT · Astro_CB · Purkinje · FC · Granule · MLI · Oligo.

### 4b. Couches / anatomie — fonds muets
| Rôle | Hex | RGB |
|---|---|---|
| ML | `#e9f1fb` | 233, 241, 251 |
| PCL | `#fdefd8` | 253, 239, 216 |
| GCL | `#efe9f7` | 239, 233, 247 |
| WM | `#eceef1` | 236, 238, 241 |
| Pie-mère | `#7f8791` | 127, 135, 145 |
| Bordure de bande | `#cfd4da` | 207, 212, 218 |

### 4c. Sémantique (DE)
| Rôle | Hex | RGB | 2ᵈ encodage |
|---|---|---|---|
| ↑ FC significatif | `#d62728` | 214, 39, 40 | — |
| ↑ BG significatif | `#2a78d6` | 42, 120, 214 | — |
| Contamination | `#e8a13a` | 232, 161, 58 | **marqueur ✕** |
| Non-significatif | `#b8b6b0` | 184, 182, 176 | petits points |
| Équivalence / null | `#5a5a57` | 90, 90, 87 | pointillés |

### 4d. Encre & chrome
| Rôle | Hex | RGB |
|---|---|---|
| Encre primaire (texte) | `#0b0b0b` | 11, 11, 11 |
| Encre secondaire | `#52514e` | 82, 81, 78 |
| Muet (ticks) | `#898781` | 137, 135, 129 |
| Grille (hairline) | `#e1e0d9` | 225, 224, 217 |
| Axe / baseline | `#c3c2b7` | 195, 194, 183 |
| Surface figure | `#ffffff` | 255, 255, 255 |

### 4e. Séquentiel (heatmaps) — 1 teinte bleu, clair→foncé (`seq_blue`)
`#cde2fb · #b7d3f6 · #9ec5f4 · #86b6ef · #6da7ec · #5598e7 · #3987e5 · #2a78d6 · #256abf · #1c5cab · #184f95 · #104281 · #0d366b`
2ᵉ contexte simultané → `seq_aqua` (`#d4f2e6 … #0e6f50`).

### 4f. Divergent (log2FC) — BG↔gris↔FC (`div_bwr`)
Pôle − `#2a78d6` · `#8fb4e6` · milieu `#f0efec` · `#e79596` · pôle + `#d62728`.

### 4g. Opacités (alpha) fixées
| Usage | Alpha |
|---|---|
| Masques de couches (overlay spatial) | **0,70** |
| Points de fond (jitter / strip) | **0,35** |
| Points ns (volcano) | **0,55** |
| Halo de mise en avant | **0,15** |
| Fills pleins (barres, chips) | **1,00** |

---

## 5. Épaisseurs de trait (pt) — **minimum 0,75 pt** (cf. §0)
| Élément | Épaisseur |
|---|---|
| Spine (axe left/bottom) | **1,0** |
| Tick major | **1,0** (longueur 3 pt, vers l'extérieur) |
| Grille (si activée) | **0,75** |
| Ligne de données | **2,0–2,5** |
| Contour de marque / barre | **0,75** |
| Contour de barre relief | **0,75** (obligatoire) |
| Barre d'erreur | **1,0** (capsize 3 pt) |
| Ligne de seuil (volcano) | **1,0** pointillés `#c3c2b7` |
| Ligne de repère / leader | **1,0** |
| Dendrogramme | **1,0** `#52514e` |
| Bordure de bande (schéma) | **1,2** `#cfd4da` |
| Pie-mère (schéma) | **3,0** `#7f8791` |

---

## 6. Marques (points, barres, lignes) — **minimum 4 pt** (cf. §0)
| Élément | Spéc |
|---|---|
| Marqueur nuage (défaut) | **6 pt** (`lines.markersize`) |
| Point ns / fond | **4 pt**, sans contour, `#b8b6b0` |
| Point mis en avant | **8 pt**, contour `#333` 0,75 pt |
| Contour de marqueur | 0,75 pt (`markeredgewidth`) |
| Largeur de barre | 0,8 unité de catégorie ; **écart 2 px** entre barres adjacentes |
| Extrémités de barre | droites, ancrées à la baseline |
| Jitter (strip) | largeur 0,15 · alpha 0,35 · **4 pt** |
| Médiane sur strip | tick épais noir, 2 pt |

---

## 7. Espacement & disposition
| Paramètre | Valeur |
|---|---|
| Gouttière inter-panneaux | 4 mm |
| Pad titre d'axe | 2,5 pt (`labelpad`) |
| Pad titre de panneau | 4 pt (`titlepad`) |
| Position lettre de panneau | coin sup-gauche, **hors** cadre, x≈−0,02 fig, y≈1,02 axes |
| Marge de données | 2 % |
| Alignement | bords de panneaux d'une rangée alignés ; baselines cohérentes |

---

## 8. Axes / ticks / grille / spines
- **Spines top & right : retirés.** Left & bottom : `#c3c2b7`, 0,6 pt.
- **Ticks** : vers l'extérieur, 2,5 pt de long, 0,6 pt de large, couleur `#52514e`.
- **Grille** : **désactivée** par défaut. Si nécessaire → hairline `#e1e0d9` 0,5 pt, **axe des valeurs uniquement**, derrière les données.
- Axes numériques : commencer à **0** pour les magnitudes (barres) ; échelle log uniquement si justifiée et annotée.
- `axes.autolimit_mode = round_numbers`.

## 9. Légendes
- **Sans cadre** (`frameon = False`).
- ≥ 2 séries → légende **présente** ; ≤ 4 séries → **aussi** labellisées en direct ; 1 série → **pas** de légende (le titre nomme).
- `handlelength = 1.2`, `handletextpad = 0.5`, `columnspacing = 1.0`.
- Position : coin le plus vide, ou au-dessus des panneaux sur une rangée.

## 10. Colorbars
- Épaisseur 6–8 % de l'axe (`fraction ≈ 0.046`), pad 0,02.
- Étiquette 7 pt, graduations 6,5 pt, **échelle affichée** (ex. 0,90–1,00), bornes explicites.
- Séquentiel `seq_blue` ; divergent `div_bwr` centré sur 0.

## 11. Barres d'échelle (figures spatiales)
- Segment plein `#0b0b0b`, épaisseur 2 pt, longueur ronde annotée (ex. **100 µm**), label 6,5 pt dessous.
- **Pas d'axes x/y** sur les vues spatiales : barre d'échelle **obligatoire** à la place.

## 12. Formatage des nombres & labels
- Décimales : cohérentes par colonne ; log2FC à **2 décimales**, corrélations à **2–3**, P en notation scientifique.
- **P** : `P = 2,5 × 10⁻⁸` (fr) / `P = 2.5e-8` (en) ; FDR = **P_adj** ; équivalence = **P_TOST**.
- Milliers : espace fine (fr) `3 570` / virgule (en) `3,570`.
- Unités toujours : `µm`, `log2`. Axe log2FC annoté « + = ↑ FC ».
- Nombres d'axes/tableaux en `tabular-nums`.

## 13. Conventions statistiques
- **DE** : Mann-Whitney U bilatéral + **BH-FDR**. Significatif = **P_adj < 0,05 ET |log2FC| > 0,25** (les deux).
- **n** reporté pour **chaque** groupe dans **chaque** panneau (ex. « FC n = 190 ; BG n = 3 570 »).
- **Équivalence** : TOST, **borne explicite** (ex. ±0,40 = ½ écart Bergmann–astro), reporter **P_TOST**.
- **Étoiles** (si employées) : `*` P<0,05 · `**` P<0,01 · `***` P<0,001 · `ns`. **Préférer le P exact.**
- Modules : z-score par gène, moyenne du panel, MWU.

## 14. Gabarits par type de figure
- **Schéma (Fig 1A)** : bandes §4b ; cellules en couleurs d'entité §4a (Bergmann/BG = `#2a78d6`, FC = `#d62728`, Purkinje = `#4a3aa7`) ; cote `≥ 45 µm` en `#c0392b` ; « not to scale » 7 pt italique.
- **Volcano** : ns `#b8b6b0` (4 pt, alpha 0,55) ; significatifs par sens (`div_bwr`) ou FC/BG ; contamination **✕** `#e8a13a` ; seuils pointillés `#c3c2b7` 1,0 pt ; labels gènes *italique* 8 pt sélectifs (peu de gènes labellisés).
- **Heatmap corrélation** : `seq_blue`, échelle 0,90–1,00, dendrogramme 0,6 pt `#52514e`.
- **Barplot marqueurs** : entités colorées, contour relief, écart 2 px, axe « expression moyenne (log2) ».
- **Strip Bergmann↔Astro** : jitter §6 ; médiane tick noir ; pôles pointillés.
- **Overlay spatial** : masques §4b alpha 0,70 ; BG points `#2a78d6`, FC étoiles `#d62728`, Purkinje `#4a3aa7` ; barre d'échelle §11.

## 15. Nomenclature & abréviations
- **Gènes** : symbole souris, *italique*, 1ʳᵉ lettre majuscule (*Gfap*, *Calb1*).
- **Groupes** : **FC** = Fañanas cell (displaced Bergmann, soma in ML) ; **BG** = Bergmann glia (PCL). Défini au 1ᵉʳ usage de chaque figure.
- **Couches** : ML, PCL, GCL, WM (développées au 1ᵉʳ usage).
- Toujours « Fañanas » avec le ñ. Unités `µm`, `log2`.

## 16. Anti-patterns → corrections
| Figure | Problème | Correction |
|---|---|---|
| `threshold_scan.png` | **double axe Y** | scinder en 2 panneaux (1 axe chacun) |
| `identity_figure.png` (A) | colormap **magma** | `seq_blue` |
| `Fig1A_*.svg` | Bergmann **teal** | Bergmann/BG → `#2a78d6`, Purkinje → `#4a3aa7`, FC → `#d62728` |

## 17. Reproductibilité
- Nommage : `FigN<lettre>_<slug>.{pdf,svg,png}` (ex. `Fig2b_classification.pdf`).
- En-tête de chaque script : `plt.style.use(...)` + `from charte_graphique.palette import ...` ; graine fixée si aléatoire.
- Planche visuelle : [`style_reference.svg`](style_reference.svg).
- Validation palette : `node dataviz/scripts/validate_palette.js "<hex,…>" --mode light --surface "#ffffff"`.
