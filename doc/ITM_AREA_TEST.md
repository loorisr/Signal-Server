# Test de `ITM_AREA_CR` (mode 14)

Branche : `test/itm-area`

## Objectif

Tester le mode « area » du NTIA ITM, déjà présent dans
`src/models/NTIA/itm_area.cpp` mais jamais appelé. Ce mode ne demande que la
distance et le delta-h, sans profil de diffraction complet.

Ajouts :

- modèle `ITM_AREA = 14`, utilisable avec `-pm 14` ;
- calcul du delta-h avec `d1thx()` sur une fenêtre sans horizon
  (`min(15·hauteur, 10 % de la distance)` à chaque extrémité) ;
- appel de `ITM_AREA_CR` avec les mêmes paramètres radio que l'ITM existant
  (climat, N0, fréquence, polarisation, epsilon, sigma, confiance,
  fiabilité), siting criteria `RANDOM`.

## Résultats

Carte : 45.0N 5.0E (point de test générique), 100 km, 1200 ppd, txh 10 m, rxh 10 m, 145 MHz, 50 W,
4 threads.

| Version | Temps 4 threads | Gain vs original |
|---|---:|---:|
| Original | 4,09 s | 1,00× |
| Exact optimisé | 3,09 s | 1,32× |
| `ITM_AREA_CR` (`-pm 14`) | **2,10 s** | **1,95×** |

Précision contre la carte exacte originale :

| Métrique | Valeur |
|---|---:|
| Écart max | 88 dB |
| p90 | 30 dB |
| p99 | 43 dB |
| p99,9 | 63 dB |
| Pixels à ≤1 dB | 5,92 % |
| Pixels à ≤5 dB | 21,69 % |

## Variantes testées

Aucune n'améliore significativement :

- fenêtre delta-h complète (0 → distance) ;
- siting criteria `CAREFUL` ou `VERY_CAREFUL` ;
- `mdvar = 12` au lieu de `SINGLE_MESSAGE_MODE`.

Les écarts restent du même ordre (p90 ≈ 30–40 dB).

## Conclusion

`ITM_AREA_CR` est ~1,5× plus rapide que l'exact, mais beaucoup trop imprécis
sur ce terrain alpin : il ignore la diffraction par le profil et prédit une
couverture très différente dans les zones d'ombre.

Il ne permet pas d'atteindre les objectifs 1 dB ou même 5 dB. Il n'est pas
fusionné dans la branche principale.

Pistes éventuelles :

- l'utiliser seulement comme prévisualisation très grossière ;
- l'utiliser en hybride : area mode sur les zones plates (faible delta-h),
  ITM profil complet ailleurs ;
- ne pas l'utiliser pour des cartes de couverture fines.
