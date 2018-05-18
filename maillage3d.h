#ifndef MAILLAGE3D_H
#define MAILLAGE3D_H
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <limits>
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <complex>

#include"C3.h"

using namespace std;

typedef double R; //Assimillation des réels aux double
typedef complex<double> C;

typedef pair<int,R3> vertexR3;
typedef pair<int,pair<int,int>> edgeR3;
typedef pair<int,tuple<int,int,int>> faceR3;
typedef pair<int,tuple<int,int,int,int>> elementR3;

class Maillage3D {

    // membres
private :
    int n_nodes ;     /*!< Nombre de noeuds du maillage*/
    int n_elems ;     /*!< Nombre d'elements du maillage*/
    int n_partitions;

    //Map des elements du maillage


    map<int,R3>* nodes;   /*!< Map des noeuds du maillage (clef) et de leur position spatiale (value)*/
    map<int,R3>* boundary_nodes; /*!< Map des noeuds du maillage situés sur le bord de celui-ci*/

    vector<int>* partitions; /*!< Set des partitions*/
    map<int,pair<int,int>>* edges; /*!< Map des elements-segments du maillage liés à leur noeuds  */

    map<int,tuple<int,int,int>>* faces; /*!< Map des elements-face explicitement définis du maillage liés à leur noeuds  */
    map<int,tuple<int,int,int>>* faces_bord;
    map<int,int>* faces_bord_partitions; /*!< Map des elements-face (clef) présents sur les bords du domaine et de sa partition correspondante (value)
                                           Note: les elements de cette map sont liés à une unique partitions, ceux-cis étant exclusivement présents sur les
                                           bords, donc liés à une seule partition*/

    map<int,tuple<int,int,int,int>>* elements; /*!< Map des elements-tetraetriques du maillage liés à leur noeuds  */
    map<int,int>* elements_partition; /*!< Map des éléments tétraèdriques du maillage et leurs partition */

    map<int,set<int>*>* neighbors; /*!< Map des noeuds du maillage et de leurs set de voisins (noeuds) respectifs*/
    map<int,set<int>*>* neighbors_elements; /*!< Map des noeuds du maillage et de leurs set de voisins (elements) respectifs*/
    map<int,set<int>*>* neighbors_faces; /*!< Map des noeuds du maillage et de leurs set de voisins (faces explicites) respectifs*/

    map<int,set<int>*>* boundary_neighbors;
    map<int,pair<int,int>>* neighbors_interface; /*!< Map des elements 3D du maillage dont au moins un sous element est dans une interface */

    map<int,pair<int,int>>* triangles_interfaces; /*!< Map des triangles (clef) et des domaines (value)s'interfacant sur ce triangle */
    map<int,tuple<int,int,int>>* triangles_interfaces_nodes; /*!< Map des triangles (clef) et ses noeuds (value) situés sur une interface du maillage */
    // methodes
public :

    Maillage3D(ifstream& file);

     vector<int>* GetPartitionsMap(){return partitions;}
     map<int,int>* GetPartitionDataMap(){return elements_partition;}
     map<int,R3>* GetNodesMap(){return nodes;}
     map<int,R3>* GetBoundaryNodesMap(){return boundary_nodes;}
     map<int,tuple<int,int,int,int>>* GetElementsMap(){return elements;}
     map<int,set<int>*>* GetNeighboursMap(){return neighbors;}
     map<int,set<int>*>* GetNeighboursElementsMap() {return neighbors_elements;}
     map<int,set<int>*>* GetNeighboursFacesMap() {return neighbors_faces;}
     map<int,tuple<int,int,int>>* GetFacesMap() {return faces;}
     map<int,tuple<int,int,int>>* GetBoundaryFacesMap() {return faces_bord;}
     map<int,pair<int,int>>* GetTrianglesInterfaceMap() {return triangles_interfaces;}
     map<int,tuple<int,int,int>>* GetTrianglesNodesInterfaceMap() {return triangles_interfaces_nodes;}


    int GetPartitionSize(int part);
    int GetElementsSize(){return elements->size();}
    int GetFacesSize() {return faces->size();}
    int GetEdgesSize(){return edges->size();}
    int GetNodesSize(){return nodes->size();}

    map<int,tuple<int,int,int,int>>::iterator GetElementsBeginIterator(){return elements->begin();}
    map<int,tuple<int,int,int,int>>::iterator GetElementsEndIterator(){return elements->end();}

    int GetElementsFirstNodeIndex(int tetr){return get<0>(elements->at(tetr));}
    int GetElementsSecondNodeIndex(int tetr){return get<1>(elements->at(tetr));}
    int GetElementsThirdNodeIndex(int tetr){return get<2>(elements->at(tetr));}
    int GetElementsFourthNodeIndex(int tetr){return get<3>(elements->at(tetr));}

    void DisplayInterfacesElements();
    void DisplayInterfacesTriangles();
    void DisplayBoundaryElements();
    void DisplayBasicData();

    //Methode d'extraction des triangles depuis les elements : boucle double sur les elements proches et recherche de O(N²) VALEURS dans des maps
    // Opération critique pour le pré-traitement du maillage
    // Envisager de modifier la méthode avec boost et une map multi-indexée pour améliorer les perfs de la recherche
    void CreateInterfaceTriangles();

    int ExportMeshAsOBJFile(string FilePath,bool verbose=0);

    pair<bool,tuple<int,int,int>> CommonTriangle (tuple<int,int,int,int>,tuple<int,int,int,int>);

    ~Maillage3D() ;
};

#endif // MAILLAGE3D_H
