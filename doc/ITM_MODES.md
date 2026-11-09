# ITM : mode exact et mode rapide (`-fast`)

Branche : `claude/itm-performance-modes-j9wxv5` (basée sur
`perf/itm-exact-deterministic`).

## Les deux modes

| Mode | Option | Résultat |
|---|---|---|
| Exact (défaut) | aucune | **identique au bit près** au code d'origine exécuté sur 1 thread, quel que soit le nombre de threads |
| Rapide | `-fast [N]` (N=10 par défaut) | approché : horizons par enveloppe convexe, delta-h recalculé tous les N points et interpolé entre deux |

`-fast 1` ne garde que l'accélération des horizons (delta-h calculé à chaque
point) : c'est un mode quasi exact. `-dh` reste disponible et se combine avec
`-fast`.

## Résultats

Cas de test : émetteur sur un site de montagne dans les Alpes, vers
les vallées et les sommets voisins ; VHF 145 MHz, 50 W, émetteur à 20 m, récepteur à
2 m, `-dbm`, 1200 ppd (Copernicus GLO-90). Référence : code d'origine
(`master`) sur **1 thread**. Temps : meilleur de 5 exécutions, 4 threads
(Xeon 2,3 GHz virtualisé ; les rapports comptent plus que les valeurs).

### Rayon 100 km (5,1 M pixels)

| Mode | Temps | Gain | ≤1 dB | ≤3 dB | p99.9 | Contour −120 dBm |
|---|---:|---:|---:|---:|---:|---:|
| Original, 4 threads | 5,84 s | 1,00× | 86,30 % ¹ | 91,02 % | 89 dB | — |
| **Exact** | 4,11 s | 1,42× | **100 %** (0 pixel) | 100 % | 0 | 0 % |
| `-fast 1` | 2,45 s | 2,38× | 99,993 % | 99,998 % | 0 | 0,001 % |
| `-fast 5` | 2,05 s | 2,85× | 99,992 % | 99,998 % | 1 dB | 0,044 % |
| **`-fast 10`** | 1,86 s | 3,14× | 99,984 % | 99,997 % | 1 dB | 0,061 % |
| `-fast 20` | 1,76 s | 3,32× | 99,944 % | 99,992 % | 1 dB | 0,082 % |
| `-fast 50` | 1,65 s | 3,53× | 99,614 % | 99,960 % | 3 dB | 0,141 % |

Sur 1 thread : original 21,9 s, exact 15,2 s, `-fast 10` 6,0 s (3,6×).

¹ Le code d'origine en multi-thread n'est pas déterministe : le pixel est
calculé par le premier thread qui l'atteint, donc par un rayon (un profil)
différent. Sur cette carte 13,7 % des pixels diffèrent de plus de 1 dB de la
version 1 thread, jusqu'à 113 dB. L'erreur de `-fast 10` est deux ordres de
grandeur plus faible que cette variabilité intrinsèque.

### Rayon 50 km (1,28 M pixels)

| Mode | Temps | ≤1 dB | Pixels différents |
|---|---:|---:|---:|
| Original, 4 threads | 1,25 s | 91,05 % | 214 735 |
| Exact | 0,91 s | 100 % | 0 |
| `-fast 1` | 0,72 s | 99,999 % | 124 |
| `-fast 10` | 0,56 s | 99,954 % | 70 238 (quasi tous à 1 dB) |
| `-fast 20` | 0,52 s | 99,750 % | 103 535 |

Vérifications supplémentaires du mode exact (0 pixel différent) : 1, 3 et
4 threads ; plaine de Bresse à 446 MHz ; modèles NTIA (`-pm 13`) et LOS.

## Ce qui a changé

### Mode exact (bit-identique)

1. **Sommes préfixées retirées du chemin exact.** Elles produisaient des
   écarts d'arrondi (7 pixels à ±1 dB) sans gain mesurable. Elles ne servent
   plus qu'en `-fast`.
2. **Correction de l'`owner_map`** (bug de la branche de départ). Un pixel
   pouvait être attribué à un rayon qui s'arrête avant de l'atteindre
   (`distance <= max_range` dans `PlotPropPath`) : il n'était alors calculé
   par personne (17 pixels de bord à −200 dBm sur 100 km). La fusion rejoue
   maintenant la condition d'arrêt du code séquentiel.
3. **Sélection des quantiles du delta-h** (`get_two_qtiles`) : comptage sans
   branchement (vectorisable), compaction des bandes par masques AVX2 (repli
   scalaire portable), largeur de bande adaptative, tri par insertion sur les
   petites bandes. Le résultat est une statistique d'ordre exacte, donc
   inchangé.
4. **`ReadPath` sans altitude** pour la passe `owner_map`, qui n'a besoin
   que des coordonnées : supprime des accès mémoire aléatoires au DEM.

Limite : `hzns` reste le premier poste du mode exact. Sur relief alpin il y
a ~50 mises à jour d'horizon RX par pixel, chacune dépendant de la
précédente (`the1 += q/d`, division incluse) ; cette chaîne ne peut pas être
raccourcie sans changer les bits. Un branch-and-bound bit-exact a été écrit
et validé (horizons identiques sur 1 M d'appels) mais n'apportait aucun
gain, ni en montagne ni en plaine : il a été retiré.

### Mode rapide (`-fast`)

1. **Horizon TX en O(1)** : par rayon, `c_i = (z_i - z_tx)/d_i` et l'argmax
   préfixe de `c_i - qc_ref·d_i`.
2. **Horizon RX en O(log n)** : l'angle vu du récepteur s'écrit comme la
   pente vers les points `(x_i, z_i - qc·x_i²)` ; le point d'horizon est le
   point de tangence de l'enveloppe convexe supérieure du profil, construite
   incrémentalement le long du rayon (deux enveloppes : flux principal et
   anticipation du delta-h, chacun monotone).
   `qc = gme/2` dépend en théorie de la hauteur moyenne du profil de chaque
   pixel ; on choisit le point avec un `qc_ref` par rayon, puis l'angle est
   recalculé avec le `qc` exact du pixel. Mesuré : 218 pixels différents au
   lieu de 105 — négligeable.
3. **Distance d'horizon exacte** : la fenêtre du delta-h dépend de
   `dl[0]`/`dl[1]`, que `hzns` accumule en float (`d -= xi`). Ces distances
   sont reproduites au bit près : table par rayon côté TX, et côté RX saut en
   O(nombre de binades) de la récurrence flottante (dans une binade, la
   récurrence est une progression arithmétique exacte). Validé contre la
   récurrence naïve sur des millions de cas. Sans cela, 2 053 pixels
   différents au lieu de 105 en `-fast 1`.
4. **Anticipation du delta-h moins chère** : l'ancre suivante ne calcule
   plus que le delta-h (`ITM_delta_h_fast`) au lieu de l'ITM complet, et sa
   valeur est réutilisée au début du segment suivant au lieu d'être
   recalculée.
5. **Sommes préfixées** pour `zsys` et `z1sq1` (O(1) au lieu de O(n)).

### Postes restants en `-fast 10`

`lrprop`/`adiff`/`alos`/`avar` et leurs `log`/`exp` (≈35 %), chargement du
DEM (≈15 %, coût fixe), lecture des profils, delta-h aux ancres.

## Reproduire

```bash
mkdir build && cd build && cmake ../src && make -j
# exact
./signalserver -dem DEM/ -lat 45.0 -lon 5.0 -txh 20 -rxh 2 -f 145 \
  -erp 50 -dbm -rt -120 -R 100 -pm 1 -o carte
# rapide
./signalserver ... -fast 10 -o carte
```

`SS_RAW_DUMP=/chemin/carte.raw` écrit la grille brute (int32) pour comparer
deux exécutions pixel à pixel.
