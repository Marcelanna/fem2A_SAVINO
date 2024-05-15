# README.txt
Utilisation de programmes de simulations de processus : fem2A_SAVINO

### Initialisation
Pour lancer les programmes, il faut se situer dans le fichier fem2A auquel le fork du git a été associé via le cmd Ubuntu.
Une fois dans le dossier du git, il est possible de :
- avoir des explications sur le fonctionnement du programme via la commande ./build/fem2A -h
- lancer les tests via la commande ./buil/fem2A -t
- lancer les simulations via la commande ./build/fem2A -s
- avoir des explications sur la simulation testée via la commande./build/fem2A -v

### Tests
Les fonctions de bases des simulations peuvent être testées indépendamment. 
Tests possibles pour fonctions des classes ElementMapping et ShapeFunctions, fonctions assemble_elementary_matrix(), local_to_global_matrix(), apply_dirichlet_boundary_conditions(), assemble_elementary_vector(), local_to_global_vector().

Chacune des fonctions est associée à une variable booléenne (nom évocateur) dans le programme main.cpp. Pour tester une fonction, il faut donc associer la valeur "true" au booléen et entrer les paramètres désirés dans l'appel de la fonction en main.cpp. 
Les détails sur les paramètres d’entrée se trouvent dans le fichier fem.h. 
La déclaration et la définition des fonctions tests sont faites en tests.h. Elles sont appel aux fonctions utilisées pour la simulation (leur déclaration se fait en fem.h et leur déclaration en fem.cpp). 

### Simulations
Les simulations utilisent les fonctions des classes ElementMapping et ShapeFunctions et elles des éléments finis. La simulation qui résout le problème de Dirichlet pur, celle qui qui résout le problème de Dirichlet avec terme source et celle qui résout le problème du sinus bump sont associées à une variable booléenne qu’il faut mettre en "true" dans qu’elles soient respectivement lancées. Les fonctions sont définies dans le fichier simu.h.

Le maillage sur lequel est effectué la simulation est issu du dossier "data" et est renseigné dans le fichier main.cpp lors de l’appel de la fonction.

### Affichage simulations
Lors du lancement de la simulation pour le cas et le maillage voulu, des fichiers .mesh et .bb sont créés dans /data/output (modifier le nom au besoin). Le fichier solution .bb est visulaisé dans le logiciel Medit.  

