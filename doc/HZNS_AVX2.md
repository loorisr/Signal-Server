# `hzns` en AVX2 — expérience

Branche : `perf/hzns-avx2`

Cette branche teste une accélération AVX2 de `hzns()`, le calcul des horizons
ITM. Elle n'est **pas destinée à être fusionnée en l'état** : le gain de
vitesse est réel, mais les résultats ne sont pas identiques au mode scalaire.

## Principe

La boucle scalaire traite un point à la fois et met à jour `the0`/`the1` dès
qu'un point passe au-dessus de la ligne d'horizon. L'idée AVX2 :

1. traiter 8 points d'un bloc avec `the0`/`the1` figés au début du bloc ;
2. calculer les 8 `q` TX et RX en AVX2 ;
3. si aucun point n'est au-dessus de la ligne (à `eps` près), le bloc est
   sauté : la boucle scalaire n'aurait rien mis à jour ;
4. sinon, rejouer les 8 points avec la boucle scalaire exacte.

Les distances `d_tx`/`d_rx` du bloc sont calculées vectoriellement
(`base + k*xi`), mais l'état scalaire `d_tx`/`d_rx` continue d'être mis à jour
par additions successives pour rester fidèle à la boucle d'origine.

## Résultats

Carte de validation : 45.0N 5.0E (point de test générique), 100 km, 1200 ppd, ITM, txh 10 m,
rxh 10 m, 145 MHz, 50 W.

| Version | 1 thread | 4 threads |
|---|---:|---:|
| Branche `perf/itm-exact-deterministic` | ~11,1–11,3 s | 3,12 s |
| Branche `perf/hzns-avx2` | **9,32 s** | **2,60 s** |

Soit environ **1,2×** sur le calcul exact.

## Validation

Comparaison avec la branche scalaire :

- 62 823 pixels différents sur 17 280 000 ;
- écart maximal : 28 dB ;
- p99 : 1 dB ;
- p99,9 : 2 dB ;
- 4 threads vs 1 thread sur la branche AVX2 : 0 pixel différent.

Les écarts apparaissent aux frontières d'ombre : de très petites différences
d'arrondi sur `the0`/`the1` suffisent à changer l'horizon retenu. Une partie
vient du code AVX2 lui-même, une autre du codegen flottant de la boucle
scalaire de repli (FMA/contraction), qui n'est pas bit-à-bit identique à la
boucle d'origine.

## Essais et variantes

- construction scalaire des tableaux `d_tx`/`d_rx` par bloc : nettement plus
  lente (≈12,4 s), la construction des tableaux annulait le gain AVX2 ;
- calcul vectoriel des distances + mise à jour scalaire exacte : 9,32 s ;
- seuil `eps` de 0,1 à 100 : le taux de repli reste autour de 17–36 % et les
  écarts ne disparaissent pas, ce qui confirme que le problème n'est pas
  seulement le seuil de détection.

## Conclusion

L'AVX2 donne ~20 % de vitesse mais ne conserve pas les résultats exacts.
Pour le mode exact, la branche `perf/itm-exact-deterministic` reste la
référence. Cette branche peut servir de base à un mode « rapide » assumé,
mais `-fast 20` offre déjà un meilleur compromis vitesse/précision.
