# Test de l'interpolation des horizons ITM (`-fasthz N`)

Branche : `perf/fast-hz`

## Principe

`hzns()` est l'un des postes les plus coûteux du mode exact. Ce mode calcule
les horizons (`dl0`, `dl1`, `the0`, `the1`) à des ancres, puis les interpole
linéairement pour les points intermédiaires. Le reste de l'ITM (profil,
delta-h, diffraction, variabilité) est calculé à chaque pixel.

Un override thread-local permet à `hzns()` de réutiliser l'horizon interpolé
sans rescanner le profil.

## Résultats

Carte : 45.0N 5.0E (point de test générique), 100 km, 1200 ppd, ITM, txh 10 m, rxh 10 m, 145 MHz,
50 W, 4 threads.

| Pas N | Temps | max | p90 | p99 | p99.9 | ≤1 dB | ≤5 dB |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 2,62 s | 91 | 10 | 30 | 50 | 61,67 % | 80,14 % |
| 20 | 2,32 s | 103 | 14 | 35 | 57 | 48,54 % | 71,33 % |
| 50 | 2,16 s | 107 | 19 | 41 | 62 | 34,18 % | 59,32 % |
| 100 | 2,06 s | 107 | 22 | 44 | 65 | 26,80 % | 51,96 % |

## Conclusion

L'interpolation des horizons est **très imprécise** : les horizons changent
de façon discontinue quand un nouvel obstacle devient dominant. Même N=10 ne
donne que 61,7 % à ≤1 dB.

Le mode n'est pas fusionné. Il confirme que `hzns` doit être recalculé
finement pour rester précis ; les approximations sur la perte, le delta-h ou
les horizons ne tiennent pas les objectifs 1/5 dB sur un terrain alpin.
