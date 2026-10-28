# Pistes d'optimisation restantes

État après les commits `dbb37f0` et `0b497b3` : le mode exact tourne autour de
11,1 s en 1 thread et 3,16 s en 4 threads sur la carte de validation
(45.0N 5.0E (point de test générique), 100 km, 1200 ppd, ITM, txh 10 m, rxh 10 m, 145 MHz, 50 W).

Répartition approximative du temps exact en 1 thread :

- `hzns` : ~4,0 s (36 %)
- sélection des quantiles delta-h : ~3,0 s (27 %)
- resampling / fit / detrend de `d1thx` : ~1,0 s (9 %)
- `lrprop` / `avar` / `adiff` / `alos` : ~3,1 s (28 %)

## 1. Fusionner les boucles de `d1thx` (exact, risque moyen)

`d1thx` fait quatre passes sur ~245 valeurs :

1. rééchantillonnage du profil ;
2. régression `z1sq1` ;
3. détrend ;
4. comptage des bandes pour les quantiles.

On peut :

- calculer les sommes de régression pendant le rééchantillonnage ;
- calculer le résidu et faire le comptage des bandes pendant le détrend.

Deux passes sur quatre disparaissent. Gain estimé : 0,8 à 1,5 s en 1 thread,
soit ~0,3 s en 4 threads.

Test effectué : la fusion complète (fit + bandes dans le détrend) est en fait
plus lente (12,1 s contre 11,1 s en 1 thread), car elle casse la
vectorisation de la boucle de détrend. La fusion du fit seul n'apporte rien
non plus. Seule la file d'attente dynamique a été conservée : elle fait
passer le 4 threads de 3,16 s à ~3,12 s, sans changer les résultats.

## 2. `hzns` en AVX2 avec resynchronisation exacte (le plus gros, délicat)

Traiter 8 points à la fois avec `the0`/`the1` figés :

- comparer les 8 `q` et détecter s'il y a une mise à jour ;
- si aucune mise à jour, avancer de 8 points ;
- sinon, repasser les 8 points en scalaire exact.

Le point critique est la validation : la vectorisation ne doit jamais
transformer un `q` positif en négatif. Un seuil de sécurité fin
(`q > -1e-6`) et un fallback scalaire doivent suffire. Gain potentiel :
1,5 à 2× sur 4 s, soit ~20 % du temps total exact.

## 3. File d'attente dynamique de rayons (facile, 5 à 10 %)

Les 4 threads traitent aujourd'hui des secteurs angulaires statiques. La
charge varie selon le terrain. Avec `owner_map`, chaque rayon est
indépendant : on peut remplacer les secteurs par une file globale de rayons
(`atomic next_ray`). Le résultat reste identique, et l'équilibrage devrait
faire passer le 4 threads de 3,16 s à ~2,9 s.

## 4. Optimisations mémoire et compilation (facile, 3 à 8 %)

- `owner_map` en `uint16_t` au lieu d'`int32_t` : ~34 Mo au lieu de 69 Mo.
- Aplatir `dem_data`/`dem_signal` en tableaux 1D contigus.
- PGO (`-fprofile-generate` puis `-fprofile-use`) : 5 à 10 % sur du code
  très branché comme ITM. LTO n'a presque rien apporté.
- Éviter le `= {}` systématique des structures `prop_type`/`propa_type`
  par appel.

## 5. Mode tolérant : interpolation adaptative de la perte

`-fast 20` n'interpole que le delta-h. La vraie accélération serait
d'interpoler la perte ITM elle-même :

- calculer la perte tous les `k` points ;
- interpoler entre deux ancres ;
- vérifier au milieu (et aux quarts) ;
- subdiviser seulement si l'écart dépasse ~0,4 dB.

Sur les zones lisses, on calculerait 3 à 5 fois moins de points ; sur les
frontières d'ombre, le raffinement garderait la précision. Gain réaliste :
2 à 3× de plus que le mode exact, avec un contrôle d'erreur.

## 6. Niveau système

- Cache de résultat indexé par un hash des paramètres : instantané et exact
  si la même requête revient.
- Daemon avec DEM en cache : évite de recharger/décoder les tuiles
  Copernicus à chaque appel.
- Si seul le contour au seuil `-rt` est utile : carte grossière puis
  raffinement autour du seuil.

## Ordre recommandé

1. File d'attente dynamique + fusion `d1thx` : exact, risque faible,
   ~10 à 20 %.
2. Interpolation adaptative de la perte : le plus gros gain si ±1 dB est
   acceptable.
3. AVX2 `hzns` : à tenter seulement si l'on veut rester strictement exact
   et gagner encore ~20 %.
4. Cache DEM/résultat : à faire si les requêtes se répètent.
