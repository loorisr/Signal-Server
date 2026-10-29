# Test de l'échantillonnage angulaire grossier (`-rays N`)

Branche : `perf/spatial-coarse`

## Principe

`-rays N` ne calcule qu'un rayon sur N. Les pixels appartenant aux rayons
sautés sont ensuite remplis à partir des voisins déjà calculés :

- plusieurs passes de dilatation 8-voisins ;
- moyenne des voisins calculés à chaque passe ;
- les pixels hors de toute zone calculée ne sont pas remplis.

C'est une première version, volontairement simple, de l'échantillonnage
spatial adaptatif.

## Résultats

Carte : 45.0N 5.0E (point de test générique), 100 km, 1200 ppd, ITM, txh 10 m, rxh 10 m, 145 MHz,
50 W, 4 threads.

| Mode | Temps | max | p90 | p99 | p99.9 | ≤1 dB | ≤5 dB |
|---|---:|---:|---:|---:|---:|---:|---:|
| Exact | 3,09 s | 1 | 0 | 0 | 0 | 100 % | 100 % |
| `-rays 2` | 1,79 s | 72 | 3 | 13 | 26 | 79,41 % | 95,09 % |
| `-rays 3` | 1,35 s | 70 | 5 | 16 | 30 | 68,94 % | 92,20 % |
| `-rays 4` | 1,13 s | 81 | 6 | 18 | 33 | 60,55 % | 89,37 % |
| `-rays 2 -fast 20` | 1,47 s | 72 | 3 | 13 | 26 | 79,32 % | 95,08 % |
| `-rays 3 -fast 20` | 1,13 s | 70 | 5 | 16 | 30 | 68,86 % | 92,19 % |
| `-rays 4 -fast 20` | 0,96 s | 81 | 6 | 18 | 33 | 60,48 % | 89,37 % |

## Conclusion

- **5 dB / 90 %** : `-rays 3` donne 92,20 % à ≤5 dB en **1,35 s**, plus
  rapide que `-fast 400 -dh 10` (2,09 s, 97,23 %). C'est le seul objectif
  amélioré par cette approche.
- **1 dB / 90 %** : `-rays 2` ne donne que 79,4 % à ≤1 dB. Insuffisant.
- **1 dB / 99 %** et **5 dB / 99 %** : non atteints.

Le remplissage par moyenne lisse les frontières d'ombre et l'échantillonnage
angulaire introduit une erreur difficile à corriger sans recalculer les
rayons manquants. Une version « plus proche voisin » ou une interpolation
angulaire entre deux rayons pourrait améliorer les percentiles, mais pas
jusqu'à 1 dB / 99 % sans revenir à un pas proche de 1.

Le mode n'est pas fusionné ; il reste une piste pour un aperçu 5 dB.
