# Optimisation du calcul ITM sur une zone

Branche : `perf/itm-exact-deterministic`

Commits de la branche :

- `ad752e8` — *Optimize ITM area coverage and make output deterministic*
- `1efb1be` — *Parallelize owner map and speed up delta-h quantile selection*

Ce document décrit les changements apportés pour accélérer le calcul ITM
sur une zone, sans utiliser `-fast` ni `-dh`, et en conservant des résultats
identiques à l'arrondi près.

## Contexte et objectif

Le calcul de couverture ITM est dominé par le coût par pixel :
`qlrpfl()` refaisait plusieurs parcours O(n) du profil terrain à chaque point
d'un rayon (`zsys`, `hzns`, `z1sq1`, `d1thx`), et chaque pixel était calculé
par le premier thread qui le réclamait.

Objectifs :

- conserver le modèle ITM exact (pas de `-fast`, pas de `-dh`) ;
- réduire le coût par pixel ;
- rendre le résultat indépendant du nombre de threads ;
- garder une base saine pour un second mode tolérant ±1 dB.

## Environnement de mesure

Les temps ci-dessous ont été mesurés sur la machine de développement, avec
4 cœurs pinnés (`taskset -c 0-3`), et non sur le N100 cible. Les valeurs
absolues seront différentes sur N100, mais les rapports donnent un ordre de
grandeur.

Paramètres de la carte de validation :

- TX : `45.0N 5.0E (point de test générique)`
- rayon : `100 km`, résolution `1200 ppd`
- modèle : ITM (`-pm 1`)
- `txh 10 m`, `rxh 10 m`, `145 MHz`, `50 W`
- sortie dBm, seuil d'affichage `-rt -120`

Commande de référence :

```bash
./build/signalserver -dem /chemin/DEM/ \
  -lat 45.0 -lon 5.0 \
  -txh 10 -rxh 10 -f 145 -erp 50 \
  -dbm -rt -120 -R 100 -pm 1 \
  -geotiff -number_threads 4 -o carte
```

## 1. Correction d'un crash en Release

Dans `d1thx()`, `d1thx2()` et `d1thx_no_interpolation()`, l'allocation était
écrite dans un `assert(...)` :

```cpp
assert((s = new float[n + 2]) != 0);
```

En Release, `NDEBUG` supprime l'assertion, donc l'allocation n'était jamais
faite et `s` restait non initialisé : segfault dans `d1thx`. L'allocation a
été sortie de l'assertion.

## 2. Sommes préfixées par rayon (`ITM_ray_ctx`)

Le profil ITM d'un point y est le préfixe `elev[2..y+1]` du même rayon.
Deux tableaux de sommes préfixées sont construits une fois par rayon :

- `p0[k] = sum(elev[2..k])`
- `p1[k] = sum(i * elev[i], i=2..k)`

Ils permettent de calculer en O(1) :

- `zsys`, la hauteur moyenne du profil, au lieu d'un parcours O(n) ;
- les sommes de fenêtre de `z1sq1()` (régression linéaire) au lieu d'un
  nouveau parcours à chaque appel.

Les accumulations se font en `double`, puis sont converties en `float` comme
avant. L'écart introduit est de l'ordre de l'arrondi flottant.

## 3. Sélection exacte des deux quantiles de delta-h

`get_two_qtiles()` doit retourner les 10e et 90e percentiles d'environ 245
valeurs. L'ancienne version faisait deux `std::nth_element` complets.

La nouvelle version :

1. calcule min/max ;
2. répartit les valeurs dans 64 bins ;
3. ne conserve que les deux bins contenant les rangs demandés ;
4. fait une sélection exacte dans ces deux bins.

Le résultat est identique à celui de `nth_element` (mêmes statistiques
d'ordre), mais avec beaucoup moins de travail par appel. Le tableau d'entrée
n'est plus modifié.

## 4. Suppression des allocations par point

`d1thx()` est appelé une fois par pixel. Il allouait un petit tableau à
chaque appel. Un buffer `thread_local` réutilisable a remplacé les
`new[]/delete[]`.

## 5. Conversion dBm / field strength

L'ancien code faisait, pour chaque pixel :

```cpp
rxp = erp / pow(10, (loss - 2.14) / 10);
dBm = 30 + 10 * log10(rxp);
```

Cette expression est affine en `loss`. Les constantes sont maintenant
précalculées une fois par rayon :

```cpp
dBm = (30 + 10*log10(erp) + 2.14) - loss;
```

Même chose pour la force de champ. Cela supprime un `pow` et un `log10`
par pixel.

## 6. Multithreading déterministe (`owner_map`)

Avant, chaque pixel était réclamé atomiquement par le premier thread qui le
rencontrait. Le résultat dépendait donc de l'ordonnancement, et une carte
4 threads pouvait différer d'une carte 1 thread.

La nouvelle approche :

1. construit une séquence globale de rayons, identique quel que soit le
   nombre de threads ;
2. parcourt les rayons dans l'ordre global pour assigner chaque pixel au
   premier rayon qui l'atteint (`owner_map`) ;
3. pendant le calcul, un rayon ne traite que les pixels dont il est
   propriétaire ;
4. une fois un pixel réclamé, son entrée est invalidée pour que le même
   rayon ne l'écrive pas une seconde fois.

Résultat : la carte est identique avec 1, 4 ou N threads.

Détail important : la séquence d'angles utilise une accumulation `volatile`.
Sans cela, `-ffast-math` peut transformer la récurrence `rad += step` en
`i * step`, ce qui arrondit différemment et change les pixels aux frontières
d'ombre.

## 7. Outil de validation `SS_RAW_DUMP`

Le GeoTIFF est colorisé en 8 bits, ce qui ne permet pas de comparer les
valeurs réelles. La variable d'environnement `SS_RAW_DUMP` écrit la grille
`dem_signal` brute en int32 :

```bash
SS_RAW_DUMP=/tmp/carte.raw ./build/signalserver ... -o carte
```

Cela a servi à vérifier pixel à pixel que les optimisations ne modifiaient
pas les résultats.

## 8. Travaux complémentaires (commit `1efb1be`)

### Parallélisation de `owner_map`

Le calcul des pixels traversés par chaque rayon est maintenant fait en
parallèle par tronçons de rayons. Les listes de pixels sont ensuite fusionnées
séquentiellement dans l'ordre global des rayons, ce qui conserve exactement
la priorité « premier rayon » et donc le déterminisme.

Gain observé sur la carte de validation : environ 0,3 s en 4 threads
(3,43 s → 3,16 s) et 0,5 s en 1 thread.

### Réutilisation des quantiles précédents

`get_two_qtiles()` mémorise les deux quantiles du point précédent. La plupart
des appels (≈98 %) ne sélectionnent plus que dans deux bandes étroites
(±8 % de l'écart interquantile) autour de ces valeurs, avec un repli exact sur
le sélecteur par histogramme si l'estimation est périmée.

Le résultat reste identique : 0 pixel différent par rapport au commit
précédent sur la carte de validation. Le gain est d'environ 0,5 s en 1 thread.

### `hzns`

Deux tentatives ont été mesurées puis abandonnées :

- traitement par blocs avec un horizon figé et repli scalaire : plus lent
  (12,3 s contre 11,1 s en 1 thread) ;
- formulation par maximum avec table d'inverses : plus lente et modifiait
  certains résultats aux frontières d'ombre.

`hzns` reste donc en version scalaire exacte.

## Résultats

| Version | 1 thread | 4 threads |
|---|---:|---:|
| Code original | 15,04 s | 4,12 s |
| Commit `ad752e8` (exact) | 11,72 s | 3,43 s |
| Commit `1efb1be` (exact) | ~11,1–11,3 s | 3,16 s |

Validation :

- **1 thread optimisé vs original** : 43 pixels différents de ±1 sur
  17 280 000, soit uniquement de l'arrondi flottant.
- **4 threads optimisé vs 1 thread optimisé** : 0 pixel différent.
- Les 43 écarts viennent des sommes préfixées et de la nouvelle sélection
  des quantiles, pas d'un changement de modèle.

## Mode tolérant ±1 dB

Le mode exact reste le défaut. Pour un calcul plus rapide avec une tolérance
d'environ 1 dB, `-fast 20` a été mesuré :

- 4 threads : environ 2,46 s sur la carte de validation ;
- 99,9 % des pixels à ≤1 dB ;
- 0,098 % des pixels à >1 dB ;
- 0,024 % à >2 dB ;
- maximum observé : 9 dB sur quelques limites d'ombre isolées.

Exemple :

```bash
./build/signalserver ... -fast 20 -o carte
```

## Limites et suites possibles

- `hzns()` reste le premier poste (~35–40 % du temps exact). Les versions par
  blocs et par maximum testées étaient plus lentes ou modifiaient les
  résultats ; une vectorisation exacte demanderait un travail dédié.
- La sélection des quantiles reste un poste significatif (~20–25 %). Les
  bandes autour des quantiles précédents ont déjà réduit le coût ; une
  vectorisation de la passe de comptage pourrait encore aider.
- `owner_map` ajoute de la mémoire temporaire : environ 36 Mo de listes de
  pixels pour 100 km en 1200 ppd, plus la carte d'occupation (~69 Mo).
  C'est acceptable sur 16 Go, mais à surveiller pour de très grandes zones
  ou du HD.
- `-fast 20` est un compromis, pas une garantie mathématique : il faut le
  valider sur les terrains réellement utilisés.
