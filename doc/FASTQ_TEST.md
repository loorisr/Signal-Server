# Test de l'interpolation quadratique du delta-h (`-fastq`)

Branche : `perf/fast-quadratic`

## Principe

`-fast N` interpole le delta-h linéairement entre deux ancres. `-fastq N`
utilise trois ancres (précédente, courante, suivante) et une interpolation
quadratique, avec clamp entre les valeurs des trois ancres pour éviter les
overshoots.

L'objectif est de permettre un pas N plus grand à précision égale.

## Résultats

Carte : 45.0N 5.0E (point de test générique), 100 km, 1200 ppd, ITM, txh 10 m, rxh 10 m, 145 MHz,
50 W, 4 threads. Référence : ITM exact original.

| Mode | N | Temps | max | p99 | p99.9 | ≤1 dB | ≤5 dB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `-fast` | 20 | 2,46 s | 9 | 1 | 1 | 99,90 % | 100 % |
| `-fastq` | 20 | 2,40 s | 9 | 1 | 1 | 99,90 % | 100 % |
| `-fast` | 50 | 2,24 s | 11 | 1 | 4 | 99,34 % | 99,98 % |
| `-fastq` | 50 | 2,20 s | 13 | 1 | 3 | 99,32 % | 99,98 % |
| `-fast` | 100 | 2,19 s | 18 | 2 | 6 | 98,13 % | 99,86 % |
| `-fastq` | 100 | 2,16 s | 18 | 2 | 6 | 98,01 % | 99,85 % |
| `-fast` | 200 | 2,17 s | 19 | 3 | 7 | 95,90 % | 99,79 % |
| `-fastq` | 200 | 2,13 s | 21 | 3 | 7 | 95,86 % | 99,75 % |
| `-fast` | 300 | 2,19 s | 23 | 3 | 6 | 94,85 % | 99,82 % |
| `-fastq` | 300 | 2,13 s | 23 | 3 | 6 | 94,75 % | 99,82 % |

## Conclusion

L'interpolation quadratique n'apporte **aucun gain significatif** : les temps
et les précisions sont dans le bruit de mesure. À N=50, le p99,9 passe de 4 à
3 dB, mais le maximum se dégrade (13 au lieu de 11 dB).

Le delta-h n'est pas la source principale d'erreur aux frontières d'ombre ;
une interpolation plus riche ne permet pas d'augmenter franchement N. Le mode
n'est pas fusionné.
