# `-fast` ou `-dh` : lequel utiliser ?

**`-fast` est nettement meilleure** : plus rapide et plus précise que `-dh`.
Recommandation : `-fast 10` par défaut, `-fast 20` pour un peu plus de
vitesse, et ne plus utiliser `-dh`.

## Mesure

Émetteur en montagne (Alpes), 145 MHz, 50 W, 20 m, rayon 50 km,
1200 ppd, 4 threads, avec `-demcache`. Référence : `master` sur 1 thread.
Temps : meilleur de 5 exécutions.

| Option | Temps | Pixels à ≤1 dB | ≤3 dB | max |
|---|---:|---:|---:|---:|
| Exact (aucune option) | 0,85 s | 100 % | 100 % | 0 |
| `-dh 100` | 1,01 s | 99,88 % | 99,998 % | 28 dB |
| `-dh 50` | 0,85 s | 98,82 % | 99,993 % | 28 dB |
| `-dh 20` | 0,73 s | 88,85 % | 99,22 % | 29 dB |
| **`-fast 10`** | **0,51 s** | **99,95 %** | 99,995 % | 15 dB |
| `-fast 20` | 0,45 s | 99,75 % | 99,97 % | 15 dB |
| `-fast 10 -dh 50` | 0,47 s | 99,10 % | 99,99 % | 28 dB |
| `-fast 20 -dh 50` | 0,45 s | 98,93 % | 99,96 % | 28 dB |

## Lecture

- **`-dh` seul n'a plus d'intérêt.** Depuis l'accélération du delta-h exact
  (sélection des quantiles), `-dh 50` n'est pas plus rapide que le mode
  exact et `-dh 100` est plus lent. Seul `-dh 20` gagne un peu de temps,
  avec 11 % des pixels à plus de 1 dB.
- **Combiner les deux n'apporte presque rien.** `-fast 10 -dh 50` gagne
  0,04 s sur `-fast 10` en dégradant la précision ; `-fast 20` seul est plus
  rapide et plus précis.

## Pourquoi

`-dh N` rend chaque calcul du delta-h plus grossier (N points
rééchantillonnés au lieu de ~245). `-fast N` en fait N fois moins et
interpole entre deux ancres, ce qui marche bien parce que le delta-h varie
lentement le long d'un rayon. `-fast` accélère aussi la recherche des
horizons, que `-dh` ne touche pas.
