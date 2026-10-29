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
| ITM area `-pm 14` | `test/itm-area` / `1f6a902` | 2,10 s | p90 30 dB, 5,9 % ≤1 dB | écarté |
| `-fastloss` adaptatif | `perf/fastloss` / `6d6ef6c` | 3,6–4,4 s | 92–96 % ≤1 dB | écarté |
| `-fastq` quadratique | `perf/fast-quadratic` / `c1c2ac3` | 2,13–2,40 s | identique à `-fast` | sans gain |
| Échantillonnage angulaire `-rays 3` | `perf/spatial-coarse` / `2ec9f81` | 1,35 s | 68,9 % ≤1 dB, 92,2 % ≤5 dB | utile pour 5 dB / 90 % |
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
