# Test du mode `-fastloss` (interpolation de la perte complète)

Branche : `perf/fastloss`

## Principe

Au lieu d'interpoler seulement le delta-h comme `-fast`, ce mode calcule la
perte ITM complète à des ancres, puis interpole la perte pour les pixels
intermédiaires. Il devrait donc supprimer `hzns`, `lrprop` et `avar` pour ces
points.

Deux variantes ont été testées :

1. **stride fixe** : perte calculée tous les N points, interpolation linéaire ;
2. **adaptatif** : un point milieu est calculé pour chaque segment ; si
   l'écart avec l'interpolation linéaire dépasse `tol`, le segment est
   subdivisé. Le mode ITM (LOS / diffraction / troposcatter) doit rester
   identique sur le segment.

Carte de validation : 45.0N 5.0E (point de test générique), 100 km, 1200 ppd, ITM, txh 10 m,
rxh 10 m, 145 MHz, 50 W, 4 threads.

## Stride fixe

| Stride | Temps | ≤1 dB | ≤5 dB | p99 | max |
|---:|---:|---:|---:|---:|---:|
| 2 | 3,04 s | 85,95 % | 95,70 % | 14 | 70 |
| 4 | 1,81 s | 73,81 % | 91,22 % | 19 | 78 |
| 8 | 1,18 s | 60,96 % | 84,80 % | 25 | 79 |
| 16 | 0,87 s | 48,06 % | 75,97 % | 31 | 86 |
| 32 | 0,71 s | 36,72 % | 65,75 % | 36 | 88 |

## Adaptatif

| Stride / tol | Temps | ≤1 dB | ≤5 dB | p99 | p99.9 | max |
|---|---:|---:|---:|---:|---:|---:|
| 8 / 0,25 | 4,38 s | 96,65 % | 99,15 % | 5 | 18 | 80 |
| 8 / 0,50 | 3,95 s | 95,17 % | 98,84 % | 6 | 20 | 80 |
| 16 / 0,25 | 4,19 s | 94,92 % | 98,65 % | 7 | 22 | 80 |
| 16 / 0,50 | 3,65 s | 92,66 % | 98,14 % | 8 | 24 | 80 |
| 32 / 1,0 | 2,90 s | 85,39 % | 95,95 % | 13 | 30 | 80 |
| 32 / 2,0 | 2,30 s | 78,61 % | 94,03 % | 16 | 32 | 83 |
| 32 / 5,0 | 1,61 s | 67,07 % | 89,60 % | 20 | 37 | 83 |

## Conclusion

Aucune configuration n'atteint **1 dB / 99 %** plus vite que l'exact :
- le stride fixe est rapide mais très imprécis ;
- l'adaptatif fin (tol 0,25) atteint 96,65 % à ≤1 dB, mais prend 4,38 s,
  soit plus lent que l'exact (3,09 s) ;
- les tolérances élevées sont rapides mais ne respectent pas les objectifs.

La raison est physique : la perte ITM présente des **discontinuités aux
frontières d'ombre**. Une interpolation, même adaptative, doit recalculer
l'ITM complet presque partout dans un terrain alpin pour rester à 1 dB.
`-fast` reste meilleur car il calcule l'ITM complet à chaque pixel et
n'interpole que le delta-h.

Le mode n'est pas fusionné. Il pourrait être utile sur des terrains très
lisses ou comme aperçu grossier, mais pas pour les cartes 1/5 dB.
