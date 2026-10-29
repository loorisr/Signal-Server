# Synthèse des expériences d'optimisation

Cette page résume les pistes testées au-delà du mode exact et des modes
`-fast` / `-dh`.

## Résultats

| Expérience | Branche / commit | Temps 4 threads | Précision | Verdict |
|---|---|---:|---|---|
| Exact optimisé | `perf/itm-exact-deterministic` | 3,09 s | 100 % ≤1 dB | référence |
| `-fast 50` | idem | 2,24 s | 99,34 % ≤1 dB | retenu |
| `-fast 100 -dh 50` | idem | 2,14 s | 96,56 % ≤1 dB, 99,84 % ≤5 dB | retenu |
| `-fast 400 -dh 10` | idem | 2,09 s | 97,23 % ≤5 dB | retenu |
| ITM area `-pm 14` | `test/itm-area` / `b0d9820` | 2,10 s | p90 30 dB, 5,9 % ≤1 dB | écarté |
| `-fastloss` adaptatif | `perf/fastloss` / `b5cfca3` | 3,6–4,4 s | 92–96 % ≤1 dB | écarté |
| `-fastq` quadratique | `perf/fast-quadratic` / `e3ba1a2` | 2,13–2,40 s | identique à `-fast` | sans gain |
| Échantillonnage angulaire `-rays 3` | `perf/spatial-coarse` / `d7334e8` | 1,35 s | 68,9 % ≤1 dB, 92,2 % ≤5 dB | utile pour 5 dB / 90 % |
| Remplissage angulaire par profils | test local | 3,11 s (`-rays 2`) | 83,4 % ≤1 dB, 95,8 % ≤5 dB | dominé par le remplissage simple |
| Interpolation des horizons `-fasthz 10` | `perf/fast-hz` / `d8f01b1` | 2,62 s | 61,7 % ≤1 dB, 80,1 % ≤5 dB | écarté |
| PGO | build local | 3,03 s / 2,37 s | identique | sans gain notable |

## Détail

### ITM area (`ITM_AREA_CR`)

- ~1,5× plus rapide que l'exact ;
- mais p90 = 30 dB contre l'ITM profil complet sur les Alpes ;
- seulement 5,9 % des pixels à ≤1 dB et 21,7 % à ≤5 dB ;
- les variantes de fenêtre delta-h, siting criteria et `mdvar` ne changent
  rien.

Conclusion : le mode area ignore la diffraction par le profil et n'est pas
utilisable pour les cartes 1/5 dB.

### `-fastloss` (interpolation de la perte complète)

Deux variantes :

- stride fixe : très rapide (1,2–1,8 s) mais 60–74 % à ≤1 dB ;
- adaptatif avec point milieu et comparaison de mode : 92–96 % à ≤1 dB,
  mais 3,6–4,4 s, donc plus lent que l'exact.

Conclusion : les discontinuités d'ombre imposent de recalculer l'ITM complet
presque partout. Une interpolation de perte ne peut pas battre `-fast` en
précision à temps égal sur ce terrain.

### `-fastq` (delta-h quadratique)

Interpolation à trois ancres avec clamp. Les temps et précisions sont dans le
bruit par rapport à `-fast` ; à N=50, p99,9 passe de 4 à 3 dB mais le maximum
se dégrade. Pas de gain exploitable.

### Termes de variabilité (`avar`)

Mesure du coût de `avar` en désactivant l'appel dans le chemin rapide :

- exact normal : 11,32 s en 1 thread ;
- sans `avar` : 11,20 s.

Soit ~0,12 s, environ **1 %** du temps exact. Une approximation ou une
interpolation de `avar` ne peut donc pas améliorer significativement le
temps ; elle n'a pas été implémentée pour ne pas dégrader la précision.

### Maths en float

Remplacement global des appels `exp/log/sqrt/pow/atan/sin/cos` par leurs
variantes float dans `itwom3.0.cc` :

- 1 thread : 11,28 s contre 11,32 s ;
- 4 threads : 3,09 s contre 3,09 s ;
- 47 pixels à ±1 dB au lieu de 43.

Aucun gain : les arguments float étaient déjà résolus en `expf`/`logf` par
les surcharges C++ ; les appels restants ne sont pas le goulot.

### Micro-optimisations ITM

Suppression du `= {}` de `prop_type`/`propv_type`/`propa_type` dans le
chemin rapide :

- 1 thread : 11,48 s ; 4 threads : 3,11 s ;
- 3 437 pixels à ±1 dB (champs non initialisés) ;
- aucun gain de temps.

Écarté. Les branches `if (debug)` et le prefetch n'ont pas été retenus non
plus : coût négligeable devant le calcul.

### iGPU

Le N100 dispose d'un iGPU Intel UHD (24 EU environ) avec OpenCL.
Portage possible mais peu attractif :

- ITM est très branché (horizons, diffraction) et peu vectorisable ;
- la mémoire est partagée avec le CPU (un seul canal sur N100) ;
- il faudrait transférer/maintenir les profils et les paramètres ;
- le gain serait probablement nul ou négatif pour l'ITM complet.

L'iGPU pourrait éventuellement servir aux modèles simples (LOS, FSPL, Hata)
ou au rendu, pas au cœur ITM. Le portage complet n'est pas recommandé.

### PGO

Profil généré sur deux cartes 50 km (exact et `-fast 20`), puis rebuild
`-fprofile-use`. Gain observé : ~2 % (3,05 → 3,03 s) et ~1 % en `-fast 20`.
Pas assez pour justifier un build PGO.

### Échantillonnage angulaire (`-rays N`)

Calcul d'un rayon sur N, puis remplissage des pixels manquants par moyenne
des voisins calculés.

- `-rays 2` : 1,79 s, 79,4 % ≤1 dB, 95,1 % ≤5 dB ;
- `-rays 3` : 1,35 s, 68,9 % ≤1 dB, **92,2 % ≤5 dB** ;
- `-rays 4` : 1,13 s, 60,6 % ≤1 dB, 89,4 % ≤5 dB.

C'est le meilleur temps observé pour **5 dB / 90 %** (1,35 s contre 2,09 s
pour `-fast 400 -dh 10`), mais le remplissage lisse les frontières d'ombre et
ne permet pas d'atteindre 1 dB / 90 % ou 5 dB / 99 %.

### Remplissage angulaire par profils

Calcul des rayons actifs en profil complet, puis interpolation entre les deux
rayons encadrants à distance égale, au lieu de la moyenne 8-voisins.

- `-rays 2` : 3,11 s, 83,4 % ≤1 dB, 95,8 % ≤5 dB ;
- `-rays 3` : 2,35 s, 76,4 % ≤1 dB, 93,9 % ≤5 dB ;
- `-rays 4` : 1,94 s, 72,0 % ≤1 dB, 92,5 % ≤5 dB.

Meilleure précision que la moyenne simple, mais le coût des profils complets
et du remplissage annule le gain. Dominé par `-rays 3` moyenne (1,35 s) pour
5 dB / 90 %.

### Interpolation des horizons (`-fasthz`)

Calcul de `dl0`, `dl1`, `the0`, `the1` aux ancres et interpolation pour les
points intermédiaires, le reste de l'ITM étant calculé à chaque pixel.

Même à N=10 : 2,62 s, 61,7 % ≤1 dB seulement. Les horizons changent de façon
discontinue quand un obstacle devient dominant ; l'interpolation est trop
grossière.

### Erreur de contour (seuil `-rt`)

Pourcentage de pixels dont le côté du seuil diffère de la référence exacte :

| Mode | Temps | seuil -120 dBm | -110 dBm | -100 dBm |
|---|---:|---:|---:|---:|
| `-fast 20` | 2,46 s | 0,08 % | 0,09 % | 0,04 % |
| `-fast 50` | 2,24 s | 0,14 % | 0,17 % | 0,13 % |
| `-fast 100 -dh 50` | 2,14 s | 0,38 % | 0,51 % | 0,33 % |
| `-rays 3` | 1,35 s | 2,07 % | 1,40 % | 0,58 % |

Si l'usage est une carte de **contour** et non un raster dBm complet,
`-rays 3` donne une frontière à ~1–2 % près en 1,35 s. `-fast 20` est presque
exact sur le contour (0,08 %) en 2,46 s.

## Ce qui reste

- **Échantillonnage spatial adaptatif** : calculer une grille polaire
  grossière, interpoler en 2D, raffiner les gradients. Seul levier
  algorithmique restant pour descendre sous 2 s, mais risque d'artefacts aux
  frontières d'ombre.
- **Cache DEM / cache de résultat** : gain important si les mêmes requêtes
  reviennent ; exact et indépendant du calcul.
- **Amélioration de `-fast` par stride adaptatif** : gain limité, car le
  plancher (~2,1 s) est dû à `hzns`/`lrprop` calculés à chaque pixel.
- **Carte de contour uniquement** : si seul le seuil `-rt` compte, une passe
  grossière puis un raffinement autour du seuil peut être très rentable.
