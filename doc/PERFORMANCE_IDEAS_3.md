# Pistes de vitesse après les modes exact / `-fast` (3ᵉ série)

Suite de `doc/ITM_MODES.md`. État de départ : mode exact bit-identique à
l'original, `-fast 10` ~3× plus rapide que l'original.

## Où part le temps (`-fast 10`, 1 thread, 50 km, site alpin)

- fin de la chaîne ITM (`lrprop`, `adiff`, `alos`, `avar` et leurs
  `exp`/`log`/`cbrt`) : ~35 % ;
- chargement du MNT (décompression GeoTIFF) : ~15 %, coût fixe ;
- lecture des profils, delta-h aux ancres, allocations : ~15 %.

## Les idées

| # | Idée | Perte de précision | Gain estimé | Statut |
|---|---|---|---|---|
| 1 | Cache des tuiles MNT décodées (`-demcache`) | aucune | ~15 % sur 50 km | **implémenté, retenu** |
| 2 | Interpolation adaptative de la perte, un point sur deux le long du rayon | réglable | 20–35 % espérés | **testé, écarté** |
| 3 | Approximations rapides de `exp`/`log`/`pow`/`cbrt` dans la chaîne ITM | ~1e-6 relatif | 5–15 % espérés | **testé, écarté** |
| 4 | Delta-h moins cher aux ancres (`-fast N -dh 50`) | faible | quelques % | déjà disponible |
| 5 | Réutiliser les buffers (`owner_map`, listes de pixels), types compacts | aucune | quelques %, plus en HD | à faire |
| 6 | Mode « contour seulement » : passe grossière puis raffinement autour de `-rt` | forte hors du seuil | ≥2× | à discuter (autre produit) |
| 7 | Processus serveur persistant gardant le MNT en mémoire | aucune | supprime le chargement | à discuter (architecture de l'API) |

Déjà écartés par les essais précédents : ITM « area », sous-échantillonnage
angulaire (`-rays`), PGO (~2 %), iGPU.

## Résultats

Mesures : 4 threads, meilleur de 5 exécutions, référence = `master` sur
1 thread. Émetteur en montagne (Alpes), 145 MHz, 50 W, 20 m.

### 1. Cache des tuiles décodées — retenu

`-demcache DIR` : à la première lecture, chaque tuile décodée par GDAL est
écrite dans `DIR` en float32 brut (petit en-tête + `ppd²` valeurs,
5,8 Mo par tuile à 1200 ppd, 52 Mo à 3600 ppd). Ensuite, le fichier brut est
lu à la place du GeoTIFF. Écriture atomique (fichier temporaire puis
`rename`), donc sûre avec des exécutions concurrentes. Le cache est ignoré si
la tuile source est plus récente.

Les altitudes sont exactement celles décodées par GDAL : sortie identique
(0 pixel différent, GeoTIFF identique octet pour octet).

| Cas | Sans cache | Avec cache | Gain |
|---|---:|---:|---:|
| Exact, 50 km | 0,97 s | 0,84 s | −13 % |
| `-fast 10`, 50 km | 0,60 s | 0,51 s | −15 % |
| Exact, 100 km | 4,15 s | 3,94 s | −5 % |
| `-fast 10`, 100 km | 1,92 s | 1,75 s | −9 % |

Le gain est fixe (~0,1–0,2 s par exécution) : il compte surtout pour les
petites zones et les requêtes répétées de l'API.

L'idée de ne lire que la fenêtre utile des tuiles a été écartée : le GeoTIFF
de sortie dessine le relief sur l'emprise de toutes les tuiles chargées, donc
l'image changerait.

### 2. Interpolation adaptative de la perte — écartée

Implémentation testée (option `-fastloss T`) : quand trois échantillons
consécutifs y, y+1, y+2 appartiennent au rayon, on calcule y et y+2 ; si les
pertes diffèrent de moins de T dB et que le mode de propagation est le même,
y+1 est interpolé au lieu d'être calculé.

| Cas | Temps | ≤1 dB | max |
|---|---:|---:|---:|
| Alpes 50 km, `-fast 10` | 0,554 s | 99,87 % | 15 dB |
| `-fast 10 -fastloss 0.3` | 0,557 s | 99,44 % | 64 dB |
| `-fast 10 -fastloss 0.5` | 0,533 s | 99,21 % | 64 dB |
| `-fast 10 -fastloss 1` | 0,514 s | 98,69 % | 64 dB |
| Bresse 446 MHz, `-fast 10` | 0,367 s | 99,88 % | 10 dB |
| `-fast 10 -fastloss 0.5` | 0,355 s | 99,40 % | 87 dB |

Gain de 3 à 7 % seulement, pour une précision nettement dégradée (erreurs
de plusieurs dizaines de dB en bord d'ombre). La perte varie trop d'un
échantillon au suivant, même en plaine, pour que l'interpolation se
déclenche souvent ; quand elle se déclenche à tort, elle efface une
transition d'ombre. Même conclusion que `perf/fastloss`. Code retiré.

### 3. Approximations mathématiques — écartée

Implémentation testée : seconde compilation de `itwom3.0.cc` dans un
namespace `itmfast` (utilisée seulement par `-fast`), avec `exp`, `log`,
`log10`, `pow` et `cbrt` remplacés par des polynômes float inline (erreur
relative ≤ 4e-6). Le mode exact restait intact.

| Cas (50 km) | libm standard | approximations |
|---|---:|---:|
| `-fast 1` | 0,636 s | 0,638 s |
| `-fast 10` | 0,503 s | 0,505 s |
| `-fast 20` | 0,435 s | 0,466 s |
| `-fast 10`, 1 thread | 1,511 s | 1,575 s |

La part de libm passe bien de 26 % à 5 % du temps, mais `adiff` et `alos`
grossissent d'autant : les fonctions de la glibc récente (variantes FMA à
tables) sont déjà aussi rapides que des polynômes inline. Aucun gain, code
retiré.

## Conclusion

Les deux pistes « approchées » ne paient pas : le coût restant est le calcul
ITM lui-même, pixel par pixel, et il ne se laisse ni interpoler ni
approximer à bon compte. Les gains qui restent sont hors du modèle :

- `-demcache` (fait) ;
- garder le MNT en mémoire dans un processus serveur (piste 7) ;
- réduire le nombre de pixels réellement calculés quand le besoin le permet
  (piste 6, contour seulement).

## Pistes suivantes (non implémentées)

### A. Sauter les pixels largement sous le seuil `-rt`

Les pixels sous le seuil ne sont pas dessinés (on voit le relief). On
calcule un point sur deux le long du rayon ; si les deux voisins sont
sous `rt − marge`, celui du milieu n'est pas calculé.

Potentiel mesuré sur la carte de référence 100 km (voisins horizontaux de
l'image, marge 10 dB) :

| Seuil `-rt` | Pixels dont les deux voisins sont sous seuil − 10 dB | Dont en réalité ≥ seuil |
|---|---:|---:|
| −120 dBm | 51,5 % | 0,05 % |
| −110 dBm | 66,5 % | 0,04 % |
| −100 dBm | 81,0 % | 0,02 % |

Soit 25 à 40 % de calculs ITM en moins en sautant un point sur deux ; gain
de temps estimé 20–30 % à 100 km (estimation, pas mesure). Le raster brut
devient faux sous le seuil : à réserver au cas où seule l'image compte.

### B. Pas d'ancre du delta-h proportionnel à la distance

En `-fast 10`, le delta-h pèse encore ~18 % du temps. Loin de l'émetteur,
un point de plus change très peu un profil de 1000 points : un pas de 10
près de l'émetteur et de ~50 à 100 km devrait donner 10–15 % de gain à
précision proche. À mesurer.

### C. ITM par lots SIMD

Traiter 8 pixels à la fois (AVX2) dans `lrprop`/`adiff`/`alos`/`avar` :
2–3× sur cette partie, mais réécriture lourde d'un code très branché.

Voir aussi `doc/FAST_VS_DH.md` pour la comparaison `-fast` / `-dh`.
