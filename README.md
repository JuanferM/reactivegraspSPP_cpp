# reactivegraspSPP_cpp
ReactiveGRASP pour le Set Packing Problem (SPP)

*Commande pour compiler*:
```bash
make
```

*Commande pour lancer le programme*:
```bash
./bin/DM2
```

*Commande pour effacer les fichiers générer lors de la compilation*:
```bash
make clean
```

**ATTENTION:** le répertoire `instances` doit se trouver dans le même répertoire
d'où le programme est lancé. Si le programme est lancé depuis `bin` alors `instances`
doit être dans `bin`. Pour utiliser le makefile les répertoires `include`, `src` et
`lib` doivent être présents. Toutes les commandes présentées dans ce README sont
supposées être lancé depuis le répértoire parent de `src` (sauf indication contraire).

### Dépendances et standard C++ utilisé
Le standard C++ utilisé dans ce projet est le C++20, assurez-vous d'avoir une version
de g++ supportant ce standard.

Il est possible de résoudre les instances avec GLPK, assurez-vous d'avoir GLPK
installé sur votre machine (l'installation de la bibliothèque **libglpk-dev** est
fortement recommandée et généralement suffisante pour utiliser GLPK dans notre cas).

Pour paralléliser GRASP, nous utilisons OpenMP, assurez-vous d'avoir OpenMP installé
sur votre machine (l'installation de **libomp-dev** est fortement recommandée et
généralement suffisante pour utiliser OpenMP dans notre cas).

Pour générer les plots, nous utilisons une version modifiée de la bibliothèque Matplot++
(https://alandefreitas.github.io/matplotplusplus/) qui a besoin des dépendances
suivantes :
* **gnuplot** (version 5.2.6+)
* **libx11-dev** (recommandé si vous utilisez le protocole de
  fenêtrage X11)

La bibliothèque Matplot++ utilise aussi d'autres dépendances mais la bibliothèque dynamique
_libmatplot.so_ (répertoire `lib`) comporte tous le nécessaire pour compiler Matplot++.
Cette bibliothèque fonctionne pour l'architecture **x86**. Un
script (`configure.sh`, utilisant la commande _cp_) a été créé pour
remplacer _libmatplot.so_ par une version supportant l'architecture **aarch** (ARM) si besoin.
Pour utiliser l'architecture x86 lancez depuis votre terminal :
```bash
./configure.sh x86
```
Pour utiliser l'architecture aarch :
```bash
./configure.sh aarch
```

### GLPK
Pour lancer GLPK sur les instances veuillez redéfinir USE_GLPK (constante préprocesseur dans le
fichier `src/main.cpp`) à true :
```c
#define USE_GLPK true
```
Sinon pour utiliser notre solution :
```c
#define USE_GLPK false
```

De plus, pour activer tous les messages de GLPK veuillez redéfinir VERBOSE_GLPK (constante préprocesseur dans
le fichier `src/main.cpp`) à true :
```c
#define VERBOSE_GLPK true
```
Sinon pour ne recevoir que les messages normaux :
```c
#define VERBOSE_GLPK false
```

### OpenMP
Pour activer la parallélisation veuillez redéfinir PARALLEL (constante
prépocesseur dans le fichier `src/main.cpp`) à true :
```c
#define PARALLEL true
```
Sinon pour exécuter les itérations séquentiellement :
```c
#define PARALLEL false
```

Pour préciser le nombre de threads maximum qu'OpenMP peut utiliser pour paralléliser GRASP veuillez
redéfinir MAX_THREADS (constante préprocesseur dans le fichier `src/main.cpp`) :
```c
#define MAX_THREADS <x>
```
Où _x_ est un nombre entier strictement positif (toute valeur incorrecte sera remplacée par 10).


### Paramètres du GRASP
#### Alpha
Pour modifier les valeur de alpha veuillez redéfinir ALPHA (constante préprocesseur dans le
fichier `src/main.cpp`) :
```c
#define ALPHA <list>
```
Où _list_ est une liste de flottant entre 0 et 1 séparés par des virgules et
entre accolades (Exemple : {0.0, 0.25, 0.5, 0.75, 1.0}).

#### Delta
Pour modifier le paramètre delta veuillez redéfinir DELTA (constante
préproceseseur dans le fichier `src/main.cpp`) :
```c
#define DELTA <x>
```
Où _x_ est un nombre.

#### Nombre d'itérations avant de mettre à jour les probabilités des alpha
Pour modifier le nombre d'itérations qu'il faut attendre avant de mettre
à jour les probabilités des valeurs alpha veuillez redéfinir PROBA_UPDATE (
constante préprocesseur dans le fichier `src/main.cpp`) :
```c
#define PROBA_UPDATE <x>
```
Où _x_ est un nombre entier inférieur ou égal à NUM_ITER (toute valeur incorrecte
sera corrigée).

#### Nombre d'itérations
Pour modifier le nombre d'itérations veuillez redéfinir NUM_ITER (
constante préprocesseur dans le fichier `src/main.cpp`) :
```c
#define NUM_ITER <x>
```
Où _x_ est un nombre entier strictement positif au moins égal à 2 (toute
valeur incorrecte sera remplacée par 2).

#### Nombre d'exécutions
Pour modifier le nombre d'exécutions (run) de GRASP veuillez
redéfinir NUM_RUN (constante préprocesseur dans le fichier
`src/main.cpp`) :
```c
#define NUM_RUN <x>
```
Où _x_ est un nombre entier strictement positif au moins égal à 1 (toute valeur
incorrecte sera remplacée par 1).

#### Nombre de points pour l'affichage des plots
Pour modifier le nombre de points utilisé pour afficher les runs de GRASP
veuillez redéfinir NUM_DIVISION (constante préprocesseur dans le
fichier `src/main.cpp`) :
```c
#define NUM_DIVISION <x>
```
Où x est un nombre entier strictement positif entre 1 et NUM_ITER (toute valeur
incorrecte sera remplacée par NUM_ITER).

#### (Plus profonde) descente
Pour effectuer des améliorations par recherche locale de type plus profonde descente
veuillez redéfinir DEEPSEARCH (constante préprocesseur dans le
fichier `src/main.cpp`) à true :
```c
#define DEEPSEARCH true
```
Sinon pour des descentes "normales" :
```c
#define DEEPSEARCH false
```

#### Mode intéractif
Les plots affichés sont tous intéractifs. L'actionnement de la molette de la
souris permet de zoomer sur les plots. Cependant, une fois que le programme
est terminé les plots deviennent statiques et il est impossible d'intéragir
avec ces-derniers. Un mode intéractif a été implémenté et permet d'intéragir avec les plots
tant que l'utilisateur n'appuie pas sur une touche depuis le terminal. Pour
activer le mode intéractif, veuillez redéfinir INTERACTIVE (constante
préprocesseur dans le fichier `src/main.cpp`) :
```c
#define INTERACTIVE true
```
Pour le désactiver :
```c
#define INTERACTIVE false
```


### Listes des instances qui ont été utilisées pour l'expérimentation
* didactic.dat
* pb_100rnd0100.dat
* pb_200rnd0100.dat
* pb_200rnd0300.dat
* pb_200rnd0900.dat
* pb_200rnd1500.dat
* pb_500rnd0100.dat
* pb_500rnd0700.dat
* pb_1000rnd0100.dat
* pb_1000rnd0500.dat
