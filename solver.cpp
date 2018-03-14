#include "solver.h"
#include "utils.h"

//Sparse Linear Algebra solver system for ComplexVectoriel (C3) Matrices
//
//Implemented solvers:
//Iterative Methods
// -Conjugate Gradient
// -GMRES
// -MinRES
// -ORTHODIR
//
// Direct Methods
// -Greedy LU/QR
// TODO : LU w/ decomposition tree
//
// Solve Ax=B
// A Hermitian Symetrical Sparse Matrix
// B Complex Vector
//
// A is written as a Morse Matrix (row compressed algorithm)
//
// About preconditionning the system:
// Preconditionners are not implemented as isin this module,
// since the A Matrix is considered as a const in the module
// however, the A matrix given to the solver can be a preconditionned matrix
// Moreover, we'll assume the A matrix is RIGHT - Preconditionned.
// -> NO IMPACT ON RESIDUAL MINIMIZATION
// -> CAREFUL WHILE EXTRACTING THE OUTPUT VECTOR
Solver::Solver(FEProblem* _fep, SolverMethod _meth,double _eps):A(_fep->GetA()),eps(_eps)
{
    result = 0;
    cout<<"Ajustement des vecteurs a la taille du probleme : "<<_fep->GetRHS().size()<<endl;
    initSolver(_fep->GetRHS().size());
    B=_fep->GetBC()+_fep->GetRHS();
    switch (_meth)
    {
    case GC:
        result = Solve_GC();
        break;
    case OrthoDir:
        wp=new VectorC3;
        Awp = new VectorC3;
        result = Solve_Orthodir();
        break;
    case GMRES:
        result = Solve_GMRES();
    default:
        break;
    }
}

Solver::Solver(MatSparseC3 _mat, VectorC3 _vec, SolverMethod _meth, double _eps):A(&_mat),B(_vec),eps(_eps)
{
    result=0;
    cout<<"Ajustement des vecteurs a la taille du probleme : "<<A->GetSize()<<endl;
    initSolver(A->GetSize());
    switch (_meth)
    {
    case GC:
        result = Solve_GC();
        break;
    case OrthoDir:
        wp=new VectorC3;
        Awp = new VectorC3;
        result = Solve_Orthodir();
        break;
    case GMRES:
        result = Solve_GMRES();
    default:
        break;
    }

}

//Implémentation de la méthode ORTHODIR
int Solver::Solve_Orthodir()
{
    cout<<"Starting resolution of the problem w/ the Orthodir Method"<<endl;
    int j,p;
    X=B; //We assume Xo=b
    A->MatVectMul(X,res); //res=AXo = Ab
    res-=B; // res = AXo- b
    w=res; //w=res
    A->MatVectMul(w,v); // v = Aw
    R3 v_n=v.Re_norms();
    w/=v_n;       // w=w/||v||
    wp[0].init(B.size());
    wp[0]=w;
    Awp[0].init(B.size());
    Awp[0]=v;
    Awp[0]/=v_n;  //Aw_0=v/||v||
    rho=-(res,Awp[0]);
    res+=Awp[0]*rho;
    X+=w*rho;
    cout<<"Iter 0: ||res|| ="<<res.Re_norm()<<endl;

    j=0;
    while (res.Re_norm()>eps && j<MAX_ITER)
    {
        w=Awp[j]; //w=Aw_p-1
        A->MatVectMul(w,v); //v=Aw

        // AtA - Orhogonalisation de la direction de descente par rapport aux directions précedentes
        for(p=0;p<j;p++)
        {
            C3 alpha=(v,Awp[p]); //a_p = (v.Aw_p)
            //cout<<"Alpha : "<<alpha.X_()<<endl;
            w-=wp[p]*alpha; // w = w-a_p*v_p
            v-=Awp[p]*alpha; // v = v-a_p*Aw_p
            //cout<<"v ["<<j<<"] : "<<v.X[0]<<endl;
            //cout<<"w ["<<j<<"] : "<<w.X[0]<<endl;
        }
        //cout<<"w  : "<<w.X[0]<<endl;
        //Orthonormalisation
        R3 v_n=v.Re_norms();
        w/=v_n;  // w_p+1= w /||v||
        wp[j+1].init(B.size());
        wp[j+1]=w;
        Awp[j+1].init(B.size());
        Awp[j+1]=v;
        Awp[j+1]/=v_n; // Aw_p+1= v /||v||

        rho=-(res,Awp[j+1]);
        X+=w*rho;
        res+=Awp[j+1]*rho;

        cout<<"Iter "<<j+1<<": ||res|| ="<<res.Re_norm()<<endl;
        cout<<res.Re_norm(0)<<" - "<<res.Re_norm(1)<<" - "<<res.Re_norm(2)<<endl;
        j++;
    }
    cout<<"Solution sur X : "<<endl;
    cout<<X.X[0]<<endl;
    cout<<"Solution sur Y : "<<endl;
    cout<<X.X[1]<<endl;
    cout<<"Solution sur Z : "<<endl;
    cout<<X.X[2]<<endl;
}



int Solver::Solve_GC()
{
    cout<<"Starting resolution of the problem w/ the Conjugate Gradient Method"<<endl;
    int p=0;
    //GC Initialization
    X=B; //We assume Xo=b
    A->MatVectMul(X,res); //res=AXo
    res-=B; // res= AXo-b
    cout<<"Iter 0 : ||res|| = "<<res.Re_norms().X_()<<endl;


    w=res;  // wo=res
    //GC iteration
    while (p<MAX_ITER) //Max iteration control
    {
        A->MatVectMul(w,v); //v=Aw

        rho=-(res,w);
        rho/=(v,w);
        //cout<<rho.X_()<<" - "<<rho.Y_()<<" - "<<rho.Z_();
        //cout<<"||rho|| = "<<rho.v_norm().X_()<<" - "<<rho.v_norm().Y_()<<" - "<<rho.v_norm().Z_()<<endl;

        X+=w*rho;
        res+=v*rho;
        if((res.Re_norm()/B.Re_norm())<=(eps*eps))
        {
            util::print_separator();
            cout<<"The solver converged."<<endl;
            cout<<"Solution sur X : "<<endl;
            cout<<X.X[0]<<endl;
            cout<<"Solution sur Y : "<<endl;
            cout<<X.X[1]<<endl;
            cout<<"Solution sur Z : "<<endl;
            cout<<X.X[2]<<endl;
            break;
        }
        p++;
        g=-(res,v)/(v,w);
        w*=g;
        w+=res; //w=res+g*w
        cout<<"Iter "<<p<<" : ||res|| = "<<res.Re_norms().X_()<<endl;
    }
}

int Solver::Solve_GMRES()
{
    cout<<"Starting resolution of the problem w/ the GMRES Method"<<endl;
    int p=0;
    R3 r;
    //GMRES Initialisation
    X=B; //We assume Xo=b
    A->MatVectMul(X,res); //res=AXo
    res-=B; // res= AXo-b
    r=-res.Re_norms();

        //Construction of the Arnoldi Basis Vectors
}

void Solver::initSolver(int _l)
{
    B.resize(_l);
    X.resize(_l);
    res.resize(_l);
    w.resize(_l);
    v.resize(_l);
    g=C3(0,0,0,0,0,0);
    rho=C3(0,0,0,0,0,0);
    MAX_ITER=100;
}
