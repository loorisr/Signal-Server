# Pistes vitesse / précision (complément)

Les modes actuels plafonnent autour de 2,1 s en 4 threads parce que `-fast`
et `-dh` ne suppriment qu'une partie du travail : `hzns`, `lrprop`, `avar` et
le profil restent calculés à chaque pixel. Pour descendre vraiment, il faut
agir sur le nombre d'appels ITM complets.

## 1. Nouveau mode d'interpolation de la perte (`-fastloss`)

`-fast N` interpole le delta-h, mais appelle quand même `hzns` + `lrprop`
à chaque pixel. `-fastloss` interpolerait la perte ITM complète :

- calculer la perte à des ancres le long du rayon ;
- interpoler entre deux ancres ;
- vérifier l'erreur au milieu (et aux quarts) ;
- subdiviser uniquement si l'écart dépasse un seuil.

On saute alors `d1thx`, `hzns`, `lrprop` et `avar` pour les points
intermédiaires.

Points d'implémentation :

- interpoler en log-distance plutôt qu'en distance linéaire ;
- forcer le raffinement quand le mode ITM change (LOS / diffraction /
  troposcatter) ;
- critère d'erreur ~0,4 dB pour viser 1 dB ;
- remplir les pixels intermédiaires seulement s'ils appartiennent au rayon
  (`owner_map`) ;
- appliquer le diagramme d'antenne par pixel, ou interpoler séparément la
  perte d'antenne.

Gain estimé : 2 à 3× par rapport à l'exact, soit ~1,2–1,6 s en 4 threads,
avec un seuil d'erreur réglable.

## 2. Améliorer `-fast` sans changer de mode

- Stride adaptatif : N=50–100 sur un rayon LOS lisse, N=5–10 près d'une
  crête ou d'un changement de mode.
- Interpolation quadratique du delta-h avec trois ancres : erreur O(h³) au
  lieu de O(h²).
- Ancres aux transitions physiques (changements d'horizon ou de mode).
- Interpolation en log-distance si la perte est interpolée.

Gain réaliste : 10 à 25 % à qualité égale.

## 3. Améliorer `-dh`

- N adaptatif selon la longueur de fenêtre et la rugosité du terrain.
- Estimation d'erreur : calculer delta-h avec N et 2N sur un sous-ensemble,
  puis affiner si l'écart dépasse un seuil.
- Échantillonnage stratifié plutôt qu'un pas régulier.

Le delta-h reste plafonné à ~245 points ; le gain est plus limité que
l'interpolation de perte complète.

## 4. Nouveau mode ITM « area » NTIA

`ITM_AREA_CR` existe déjà dans `src/models/NTIA/itm_area.cpp` mais n'est
jamais appelé. Il ne demande que la distance et le delta-h, sans profil
complet ni `hzns` :

- coût par pixel potentiellement O(1) ;
- précision à valider contre l'ITM profil complet ;
- peut être combiné avec `-dh 50`.

Bon candidat pour un mode 5 dB / 90–99 % très rapide ou une prévisualisation.

**Test effectué (branche `test/itm-area`)** : `-pm 14` tourne en ~2,1 s en
4 threads (1,5× plus vite que l'exact), mais l'écart p90 avec l'ITM profil
complet est de **30 dB** ; seulement 5,9 % des pixels sont à ≤1 dB et 21,7 %
à ≤5 dB. Les variantes (fenêtre complète, siting criteria, mdvar) ne changent
rien. Le mode est donc écarté pour les objectifs 1/5 dB ; il ne reste
utilisable que comme aperçu très grossier.

## 5. Échantillonnage spatial adaptatif

Calculer une grille polaire grossière (un rayon sur deux, un point sur deux),
interpoler en 2D, puis raffiner les zones à fort gradient ou proches du seuil
`-rt`. Gain potentiel 2–4×, mais risque d'artefacts aux frontières d'ombre.

## 6. Optimisations transverses

- PGO (`-fprofile-generate` puis `-fprofile-use`) : 5–10 %.
- Cache DEM / cache de résultat : gain important si les requêtes se répètent.
- `owner_map` en `uint16_t` et DEM aplati : quelques pourcents.
- AVX2 `hzns` exact : ~20 % de gain possible, mais les essais actuels ne sont
  pas bit-exacts à cause du codegen flottant.

## Recommandation

1. `-fastloss` adaptatif : meilleur gain vitesse/précision.
2. Stride adaptatif dans `-fast` : peu risqué, 10–25 %.
3. Benchmark `ITM_AREA_CR` : peut donner un mode très rapide à 5 dB.
4. PGO + cache DEM/résultat : gains faciles pour un usage répété.
