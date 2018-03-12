#ifndef SOLVER_H
#define SOLVER_H


#include "RNM/RNM.hpp"
#include "C3.h"
#include "MatSparseC3.h"
#include "feproblem.h"

enum SolverMethod {GC=1,OrthoDir=2,GMRES=3};

class Solver
{
public:
    Solver(FEProblem*, SolverMethod,double);

    Solver(MatSparseC3, VectorC3, SolverMethod, double);

    int Solve_GC();
    int Solve_Orthodir();
    int Solve_GMRES();

    void initSolver(int);
private:
    MatSparseC3* A;
    VectorC3 B;
    VectorC3* X_legacy; //previous intermediary solutions
    VectorC3 X; //current solution
    VectorC3 res; //current residual
    VectorC3 w; //current direction
    VectorC3 v;
    VectorC3* Awp; //previous directions
    double eps;
    C3 g;
    C3 rho;
};

#endif // SOLVER_H
