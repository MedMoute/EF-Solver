#ifndef FEPROBLEM_H
#define FEPROBLEM_H

#include "RNM/RNM.hpp"

#include "maillage3d.h"
#include "rhs.h"
#include "bc.h"
#include "MatSparseC3.h"
#include "muparserx/parser/mpParser.h"

#include "parametersdialog.h"


#include <QObject>

using namespace mup;


class FEProblem : public QObject
{
    Q_OBJECT
public:
    FEProblem(Maillage3D* _maillage,RHS* _rhs,BC* _bc,int _op,bool _unit);
    ~FEProblem(){delete A;}
    void ClearA(){delete A;}
    C3 CalcIntOnTriangle(Maillage3D* _maillage,BC* _bc,int idx,int tri_idx);
    C3 CalcIntOnTetrahedron(Maillage3D* _maillage,RHS* _bc,int idx,int tetra_idx);

    void BuildPoissonOperator();
    void BuildParameterPoissonOperator(double lambda);
    void BuildPartitionnedParameterPoissonOperator(std::map<int, double> lambdas,std::map<int,int>* partition_data);

    void BuildHeatOperator();
    void BuildParameterHeatOperator(double lambda);
    void BuildPartitionnedParameterHeatOperator(std::map<int,double> lambdas, std::map<int, int> *partition_data);

    void BuildHelmholtzOperator();
    void BuildParameterHelmholtzOperator(double omega);
    void BuildPartitionnedParameterHelmholtzOperator(std::map<int, double> omegas,std::map<int,int>* partition_data);

    //Méthodes de debug et de visualisation des matrices elementaires
    void DisplayElementaryMassMatrix(int n_elem);
    void DisplayElementaryStiffnessMatrix(int n_elem);

    VectorC3 GetBC(){return G;}
    VectorC3 GetRHS(){return F;}
    MatSparseC3* GetA(){return A;}
public slots:
    void ProcessParamData();
private:
    VectorC3 F;
    VectorC3 G;
    MatSparseC3* A;

    Maillage3D* maillage;

    ParametersDialog* dialog;

    int op;
    bool unitary;
};

#endif // FEPROBLEM_H
