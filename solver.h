#ifndef SOLVER_H
#define SOLVER_H


#include "RNM/RNM.hpp"
#include "C3.h"
#include "R3.h"
#include "MatSparseC3.h"
#include "feproblem.h"

enum SolverMethod {GC=1,MinRES=2,OrthoDir=3};

enum ErrorNorm {infinity=1,L1=2,L2=3};

class Solver
{
public:
    Solver(FEProblem*, SolverMethod,double);

    Solver(MatSparseC3, VectorC3, SolverMethod, double);

    void initSolver(int);

    int Solve_GC();
    int Solve_Orthodir();
    int Solve_GMRES();
    int Solve_MinRES();
    //DO NOT USE THIS METHOD FOR THE DIRECT SOLVER SINCE THE MATRIX POINTER IS UNSTABLE
    //USE THE EXPLICIT METHOD BELOW
    R3 computeError(ErrorNorm norm=ErrorNorm::infinity);

    R3 computeError(MatSparseC3 A, VectorC3 B,ErrorNorm norm=ErrorNorm::infinity);

    void displaySolution();

    int GetRes(){return result;}
    VectorC3 GetSolution(){return X;}
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

    bool verbose;//Display(or not) the final output
};

#endif // SOLVER_H
