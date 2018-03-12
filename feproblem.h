#ifndef FEPROBLEM_H
#define FEPROBLEM_H

#include "RNM/RNM.hpp"

#include "maillage3d.h"
#include "rhs.h"
#include "bc.h"
#include "MatSparseC3.h"
#include "muparserx/parser/mpParser.h"

using namespace mup;


class FEProblem
{
public:
    FEProblem(Maillage3D* _maillage,RHS* _rhs,BC* _bc);
    C3 CalcIntOnTriangle(Maillage3D* _maillage,BC* _bc,int idx,int tri_idx);
    C3 CalcIntOnTetrahedron(Maillage3D* _maillage,RHS* _bc,int idx,int tetra_idx);

    void BuildPoissonOperator();
    void BuildHeatOperator();
    void BuildHelmholtzOperator();
    //Méthodes de debug et de visualisation des matrices elementaires
    void DisplayElementaryMassMatrix(int n_elem);
    void DisplayElementaryStiffnessMatrix(int n_elem);

    VectorC3 GetBC(){return G;}
    VectorC3 GetRHS(){return F;}
    MatSparseC3* GetA(){return A;}
private:
    VectorC3 F;
    VectorC3 G;
    MatSparseC3* A;

    Maillage3D* maillage;
};

#endif // FEPROBLEM_H
