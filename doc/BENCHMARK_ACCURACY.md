# Benchmark vitesse / précision

Mesures effectuées sur la branche `perf/itm-exact-deterministic`, sur la
machine de développement avec 4 cœurs pinnés (`taskset -c 0-3`). Le N100
donnera des temps absolus différents ; les rapports sont indicatifs.

## Carte de validation

- TX : 45.0N 5.0E (point de test générique)
- rayon : 100 km, 1200 ppd
- modèle : ITM (`-pm 1`)
- txh 10 m, rxh 10 m, 145 MHz, 50 W
- sortie dBm, seuil `-rt -120`
- 5 102 169 pixels calculés

## Référence et métriques

La référence est la carte exacte du code original (`af2c672`, avec la
correction du crash Release). Pour chaque configuration :

- `≤ 1 dB` / `≤ 5 dB` : pourcentage de pixels dont l'écart absolu avec la
  référence est inférieur ou égal au seuil ;
- `p90`, `p99`, `p99.9` : percentiles de l'écart absolu ;
- `max` : écart absolu maximal observé ;
- temps : médiane de 3 exécutions en 4 threads.

## Résultats

| Solution | Options | Temps 4 threads | Gain vs original | max | p90 | p99 | p99.9 | ≤ 1 dB | ≤ 5 dB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Original | aucun | 4,09 s | 1,00× | 0 | 0 | 0 | 0 | 100 % | 100 % |
| Exact optimisé | aucun | 3,09 s | 1,32× | 1 | 0 | 0 | 0 | 100 % | 100 % |
| 1 dB / 99 % | `-fast 50` | 2,24 s | 1,83× | 11 | 1 | 1 | 4 | **99,34 %** | 99,98 % |
| 1 dB / 90 % + 5 dB / 99 % | `-fast 100 -dh 50` | 2,14 s | 1,91× | 19 | 1 | 3 | 7 | **96,56 %** | **99,84 %** |
| 5 dB / 90 % | `-fast 400 -dh 10` | 2,09 s | 1,96× | 19 | 3 | 8 | 11 | 70,87 % | **97,23 %** |
| Marge 1 dB | `-fast 20` | 2,46 s | 1,66× | 9 | 0 | 1 | 1 | **99,90 %** | 100 % |

## Lecture des seuils

Selon l'interprétation de « 1 dB / 99 % » :

- **99 % des pixels à ≤1 dB** : `-fast 50` suffit (99,34 %) ;
- **p99 ≤1 dB** (le 99e percentile des écarts) : `-fast 50` et `-fast 20`
  sont à p99 = 1 dB ; `-fast 20` a en plus p99.9 = 1 dB et max 9 dB.

Pour **5 dB / 99 %**, `-fast 100 -dh 50` donne 99,84 % de pixels à ≤5 dB.
Pour **5 dB / 90 %**, `-fast 400 -dh 10` donne 97,23 % à ≤5 dB.

## Remarques

- Les gains restent modestes (1,8 à 2,0× vs original) parce que le mode
  approximé ne supprime que le recalcul du delta-h ; le reste du calcul ITM
  (horizons, profils, variabilité) reste inchangé.
- La branche AVX2 expérimentale n'est pas incluse : elle est plus lente que
  `-fast 20` (2,60 s contre 2,46 s en 4 threads) et moins précise.
- Une interpolation adaptative de la perte ITM complète, avec contrôle
  d'erreur, reste la piste pour descendre plus bas en temps.
