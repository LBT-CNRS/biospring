# Revue technique et améliorations de BioSpring

## Périmètre

Cette revue porte sur le cœur C++ de BioSpring : réseau de ressorts, interactions non liées, recherche de voisins, sonde interactive, corps rigides, entrées/sorties NetCDF/PDB/PQR/CSV et exécution OpenMP.

Le dépôt examiné contient environ 190 fichiers C/C++ pour 34 000 lignes, avec les sous-systèmes principaux suivants :

- `spn/` : particules, ressorts et boucle de simulation ;
- `topology/` et `reduce/` : représentation de topologie et réduction gros grain ;
- `forcefield/` : potentiels ressort, électrostatique, stérique, IMPALA et hydrophobe ;
- `interactor/` : MDDriver et FreeSASA pour les interactions externes ;
- `rigidbody/` : dynamique des corps rigides ;
- `opencl/` : implémentation GPU distincte ;
- `IO/` : topologies et trajectoires.

## Corrections appliquées

### 1. Calcul des ressorts déterministe et sûr avec OpenMP

Le calcul parallèle précédent modifiait directement les forces des particules depuis plusieurs threads. Deux ressorts partageant une particule pouvaient donc écrire simultanément dans le même `Vector3f`.

Le nouveau chemin :

1. calcule en parallèle une contribution de force indépendante par ressort ;
2. accumule ensuite les forces et les énergies dans un ordre fixe ;
3. réutilise un tampon de travail entre les pas.

Cela supprime la course de données et rend forces et énergie de ressort indépendantes du nombre de threads.

### 2. Calcul des forces particulaires reproductible

Chaque thread ne modifie désormais que la particule qui lui est assignée. Les énergies sont sommées dans l’ordre des particules après la région parallèle, ce qui évite les variations dues aux réductions flottantes OpenMP.

L’agrégation des forces et couples de corps rigides est également effectuée hors de la boucle parallèle, car les accumulateurs d’un même corps sont partagés.

### 3. Sonde interactive corrigée

La sonde était auparavant :

- modifiée simultanément par plusieurs threads ;
- intégrée une fois par particule dynamique au lieu d’une fois par pas ;
- susceptible d’effacer la force de la particule courante ;
- comptée plusieurs fois dans les énergies ;
- dupliquée dans les sorties PDB/PQR.

Elle est maintenant remise à zéro une fois, alimentée séquentiellement par les interactions, intégrée exactement une fois par pas, puis synchronisée vers sa représentation dans le tableau de particules. Les interactions ordinaires ignorent explicitement cette représentation afin d’éviter le double comptage.

### 4. Recherche de voisins mise à jour pendant la dynamique

Les listes cellulaires étaient construites pendant `setup()` puis conservées alors que les particules se déplaçaient. Une particule changeant de cellule pouvait perdre ou conserver de faux voisins.

Les structures stérique, électrostatique et hydrophobe sont maintenant reconstruites avant chaque calcul public des forces particulaires. Les lectures des cellules ont aussi été rendues constantes et n’utilisent plus `unordered_map::operator[]` pendant les régions parallèles.

Cette correction privilégie la justesse. Une liste de Verlet avec peau est recommandée pour retrouver de meilleures performances à grande échelle.

### 5. Stabilité des références et pointeurs de ressorts

Les particules mémorisaient des pointeurs vers les éléments d’un `std::vector<Spring>`. Une réallocation du vecteur invalidait ces pointeurs.

Après une réallocation, les tables de voisinage élastique sont désormais reconstruites. L’ajout de particules après la création des ressorts est refusé, et une place est réservée pour la sonde avant que les ressorts ne conservent des références vers les particules.

### 6. Initialisation et remise à zéro complètes

Les champs hydrophobes, identifiants, énergies, positions précédentes et autres membres de `Particle` et `ParticleProperty` sont explicitement initialisés. L’énergie hydrophobe est remise à zéro et exposée dans le journal et la sortie CSV.

`SpringNetwork::clear()` nettoie maintenant les listes dérivées, recherches de voisins, ressorts, tampon de forces, sonde et vecteur d’insertion afin d’éviter les références résiduelles.

### 7. Activation effective du potentiel hydrophobe

La méthode `_setupHydrophobic()` existait mais n’était jamais appelée par `setup()`. Elle est maintenant intégrée à l’initialisation, avec validation du cutoff.

### 8. Corrections des entrées/sorties et de la mémoire

- Les tableaux alloués avec `new[]` sont libérés avec `delete[]`.
- Les tampons de taille nulle possèdent maintenant un état défini.
- Le tableau « nombre de ressorts par particule » est alloué même lorsqu’il n’existe aucun ressort.
- Les tampons peuvent être réinitialisés sans fuite et ne sont plus copiables accidentellement.
- La surface accessible au solvant est réellement copiée dans le tampon NetCDF.
- La lecture NetCDF conserve simultanément surface accessible et énergie de transfert IMPALA au lieu d’écraser la première propriété avec la seconde.
- Le fichier NetCDF ouvert est possédé par un `std::unique_ptr`.
- Les enregistrements factices supplémentaires de sonde ont été retirés des sorties PDB et PQR.

### 9. Tenseur d’inertie des corps rigides

Le calcul utilisait `=+` au lieu de `+=`, ce qui remplaçait les contributions au lieu de les accumuler. Plusieurs signes diagonaux et hors diagonale étaient aussi incorrects. Le tenseur suit maintenant la forme symétrique standard.

### 10. Robustesse d’API

- `setInsertionVector(aa1, aa2)` utilise réellement ses arguments.
- Les ressorts vérifient les indices et refusent une auto-liaison.
- Une liste vide de contraintes ne provoque plus de division par zéro.
- Les positions infinies sont détectées en plus des NaN, avec l’identifiant de la particule fautive.

## Validation effectuée

Les contrôles suivants ont réussi :

- vérification syntaxique C++20 des fichiers modifiés, avec et sans OpenMP selon le fichier ;
- recherche de voisins après changement de cellule ;
- lectures concurrentes d’une recherche de voisins constante ;
- construction, destruction et réinitialisation des tampons sous AddressSanitizer et UndefinedBehaviorSanitizer ;
- calcul élémentaire d’un ressort sous AddressSanitizer ;
- test intégré du cœur, compilé sans la couche NetCDF, couvrant :
  - égalité exacte des forces et de l’énergie de ressort avec 1 et 4 threads ;
  - stabilité des pointeurs après réallocations du tableau de ressorts ;
  - mise à jour du voisinage après déplacement ;
  - intégration unique de la sonde et bilan des forces ;
  - initialisation du réseau hydrophobe ;
  - système sans ressort ;
  - absence de sonde PDB dupliquée.

La chaîne CMake utilise `find_package()` et des cibles importées. NetCDF C et C++4 sont recherchés séparément : les packages config officiels sont privilégiés, puis les modules acceptent des chemins explicites vers les en-têtes et bibliothèques. OpenMP reste facultatif : il est recherché et lié via `OpenMP::OpenMP_CXX` uniquement lorsque `OPENMP_SUPPORT=ON`.

## Recommandations pour l’étape suivante

### Priorité 1 — voisinage avec peau de Verlet

Reconstruire trois grilles à chaque pas garantit la justesse mais coûte `O(N)` à chaque itération. Une liste de Verlet commune, reconstruite lorsque le déplacement maximal dépasse la moitié d’une peau configurable, réduirait fortement ce coût.

### Priorité 2 — noyaux de paires uniques

Les interactions non liées sont encore évaluées deux fois, une fois depuis chaque particule. Une liste de paires `i < j`, combinée à des tampons de forces par thread ou à une coloration du graphe, permettrait d’appliquer explicitement la troisième loi de Newton sans course de données.

### Priorité 3 — couche de données orientée calcul

`Particle` est une structure riche de type AoS. Pour les noyaux CPU SIMD et GPU, une représentation SoA séparant positions, forces, vitesses et paramètres réduirait les défauts de cache et faciliterait la vectorisation. La couche objet pourrait rester l’API de contrôle et de sérialisation.

### Priorité 4 — véritable couplage multirésolution

La réduction gros grain est présente, mais le cœur ne formalise pas encore une zone hybride tout-atome/gros-grain. Il faudrait expliciter :

- la table de projection atomes ↔ sites gros grain ;
- la conservation de quantité de mouvement lors des transferts ;
- une fonction de pondération spatiale ou topologique ;
- l’interpolation énergie/force et, si nécessaire, la force thermodynamique de compensation ;
- le traitement des interactions traversant la frontière de résolution.

### Priorité 5 — sous-espace de modes normaux interactif

Pour un réseau élastique piloté par modes normaux :

- construire une Hessienne creuse du réseau autour d’une conformation de référence ;
- extraire seulement les modes de basse fréquence avec un solveur itératif ;
- intégrer les amplitudes modales comme degrés de liberté ;
- projeter les forces interactives sur le sous-espace modal ;
- permettre la mise à jour locale ou périodique de la base lorsque la conformation s’éloigne de la référence.

### Priorité 6 — intégrateurs et invariants physiques

Euler explicite est simple mais peu stable pour des ressorts raides. Ajouter Velocity Verlet, un contrôle de pas fondé sur la fréquence maximale et des tests de dérive d’énergie améliorerait la fiabilité des simulations non dissipatives.

### Priorité 7 — unifier CPU et OpenCL

L’implémentation OpenCL est une branche de calcul distincte. Une description commune des paramètres et des noyaux, avec tests numériques CPU/GPU sur les mêmes systèmes, limiterait les divergences fonctionnelles.

## Conclusion

Le patch fourni cible d’abord les erreurs pouvant produire des trajectoires fausses, non déterministes ou des corruptions mémoire. Il constitue une base plus sûre pour une seconde étape orientée performances et couplage multirésolution formel.
