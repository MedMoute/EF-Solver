#ifndef SOLVER_H
#define SOLVER_H


#include "RNM/RNM.hpp"
#include "C3.h"
#include "MatSparseC3.h"
#include "feproblem.h"

enum SolverMethod {GC=1,OrthoDir=2,GMRES=3,MinRES=4};

class Solver
{
public:
    Solver(FEProblem*, SolverMethod,double);

    Solver(MatSparseC3, VectorC3, SolverMethod, double);

    int Solve_GC();
    int Solve_Orthodir();
    int Solve_GMRES();
    int Solve_MinRES();

    void initSolver(int);
    int GetRes(){return result;}
private:
    MatSparseC3* A;
    VectorC3 B;
    VectorC3* X_legacy; //previous intermediary solutions
    VectorC3 X; //current solution
    VectorC3 res; //current residual
    VectorC3 w; //current direction
    VectorC3 v;
    VectorC3* Awp; //previous A*directions
    VectorC3* wp; //previous directions
    VectorC3* vp;//MinRES Lanczos vectors
    double eps;
    C3 g;
    C3 rho;
    int result;
    int MAX_ITER;//nombre maximal d'iteration
    int TRUNC; // size of the truncated method (not implemented yet)
};

#endif // SOLVER_H
