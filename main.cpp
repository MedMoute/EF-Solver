/**
 * \file main.cpp
 * \brief Programme de résolution de problèmes linéaire par la méthode des Elements Finis.
 * \author Mehdi Ennaïme
 * \version 0.1
 * \date 20 Avril 2018
 *
 * Programme réalisant la résolution d'un problème d'EDPs vectorielles linéaires par le biais de la méthode des élements finis 3D.
 * Ce programme utilise des elements finis P1
 *
 */

#include "map_mainwindow.h"
#include <QApplication>

#include "C3.h"
#include "utils.h"
/*! \namespace std
*
* \brief Espace de nommage standard
*
 * Utilise le namespace de la bibliothèque standard
 */
using namespace std;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MAP_MainWindow w;
    w.show();

   return a.exec();
}


// Le contenu qui suit sert uniqument à la génération de la page d'index du documentation Doxygen
/*! \mainpage Accueil
 *
 * \section intro_sec Introduction
 *
 * L'objectif de ce projet est la réalisation d'un solveur de problème de type EDP linéaires.
 * Plus particulièrement, l'objectif principal est la résolution du problème dit du Micro-Onde, dont on fait ici le rappel de l'énoncé:
 *
 * \subsection a Enoncé du problème
 * En version Helmholtz il s’agit de résoudre un probleme du type:
 * Trouver u une fonction définie sur l’ouvert \Omega a valeur dans C tel que
 *
 *
 * ω^2µu + ∇.1/ε∇u = 0 ∈ Ω
 * avec des données aux bord imaginaires.
 *
 *
 * Le couplage avec la température θ, se fait comme suivant:
 *
 * −∇.K∇θ = f
 *
 * dans l’objet a cuire, et où f est nulle hors de l’objet a cuire, et est égal a uu dans l’objet a cuire.
 *
 * Le projet est la résolution de 2 équations découples avec les données suivantes :
 *
 *  ε = 1 dans l’air et ε = 4 dans la viande.
 *
 *
 * On pourra prendre des conditions de Dirichlet sur l’objet a cuire.
 *
 * \subsection Problématiques
 *
 * En étudiant l'énoncé ,on constate que ce projet comporte plusieurs points importants, liés à dés problématiques variées sur le plan de
 * l'analyse numérique, de l'analyse matricielle ou encore de l'informatique.
 *
 * - Le problème est tridimensionnel, cela implique une croissance cubique de la taille des matrices par rapport à la dimension du problème.
 * Cela implique que la résolution des systèmes linéaires ne peut se faire par le biais d'un algorithme de résolution directe (décompositions QR/LU),
 * on préviligiera donc des méthodes itératives.
 *
 * - Le problème est vectoriel, mais l'opérateur est unique, i.e. chaque résolution comporte 3 "seconds membres" correspondant chacuns à l'une des dimensions
 * dans laquelle le problème est résolu.
 *
 * - La matrice correspondant au problème discrétisé n'est pas nécessairment définie, elle est en revanche symétrique. Cela implique un choix préférentiel
 * pour la résolution du système : la méthode MINRes
 *
 * - Quant à la gestion des données d'entrée/sortie, l'énoncé reste vague, mais implique que le programme doit être capable de gérer les partitions de maillages.
 *
 *
 */
