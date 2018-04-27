#include "solver.h"
#include "utils.h"

//Sparse Linear Algebra solver system for ComplexVectoriel (C3) Symetrical Matrices
//
//Implemented solvers:
//Iterative Methods
// -Conjugate Gradient
// -MinRES  (NOT WORKING)
// -ORTHODIR (WORKS - UNSTABLE)
//
// Direct Methods
// TODO : LU w/ decomposition tree
//        Dense LU ? -> Matrix is stored sparsely, probably bad idea
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
        delete wp,Awp;
        break;
    case MinRES:
        wp = new VectorC3;
        vp = new VectorC3;
        result = Solve_MinRES();
        delete wp,vp;
        break;
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
        delete wp,Awp;
        break;
    case MinRES:
        wp = new VectorC3;
        vp = new VectorC3;
        result = Solve_MinRES();
        delete wp,vp;
        break;
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
    C3 v_n=v.Re_norms();
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
        for(p=0;p<=j;p++)
        {
            C3 alpha=(v,Awp[p]); //a_p = (v.Aw_p)
            //cout<<"Alpha : "<<alpha<<endl;
            w-=wp[p]*alpha; // w = w-a_p*v_p
            v-=Awp[p]*alpha; // v = v-a_p*Aw_p
            //cout<<"(w,w["<<p<<"])"<<(w,wp[p])<<endl;
            //cout<<"(v,Aw["<<p<<"])"<<(v,Awp[p])<<endl;
            //cout<<"w ["<<j<<"] : "<<w.X[0]<<endl;
        }
        //Orthonormalisation
        C3 v_n=v.Re_norms();
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
        //cout<<"Res : "<<res.Re_norms()<<endl;
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

int Solver::Solve_MinRES()
{
    cout<<"Starting resolution of the problem w/ the MinRES Method"<<endl;
    int p=0;
    C3 alpha;
    R3 norm_alpha;
    C3 r[2],t[3],c[2],s[2],a[2],b,comp;
    //MinRES Initialisation
    X=B; //We assume Xo=b
    A->MatVectMul(X,res); //res=AXo
    res-=B; // res= AXo-b
    r[0]=-res.Re_norms();
    //For the MinRES Method the direction vectors are the Lanczos Basis Vectors
    vp[0].init(B.size());
    vp[0]=res*r[0].comp_inv();
    comp=r[0];
    //Construction of the Lanczos Basis Vectors
    while (p<MAX_ITER || R(comp.norm())>=eps)
    {
        cout<<"Iteration n"<<p<<" - Residu : "<<comp<<endl;
        vp[p+1].init(B.size());
        wp[p].init(B.size());
        cout<<"init done"<<endl;
        A->MatVectMul(vp[p],w); // w=A*vp
        t[0]=t[2];
        w-=vp[p]*t[0];
        t[1]=(w,wp[p]);
        w-=vp[p]*t[1];
        t[2]=w.Re_norms();
        vp[p+1]=w*t[2].comp_inv();
        //Application of the two previous Givens Rotations
        r[0]=-t[0]*s[0];
        a[0]=t[0]*c[0];
        cout<<t[2]<<endl;
        cout<<"givens done"<<endl;
        //
        r[1]=c[1]*a[0]-s[1]*t[1];
        a[1]=s[1]*a[0]+c[1]*t[1];
        //Update c[] and s[]
        c[0]=c[1];
        s[0]=s[1];
        c[1]=C3((t[2]/a[1])+C3(1,0,1,0,1,0)).comp_inv().sqrt_();
        s[1]=-c[1]*t[2]/a[1];
        //
        r[2]=c[1]*a[1]-s[1]*t[2];
        b=c[1]*r[2];
        comp=s[1]*r[2];
        cout<<comp<<endl;
        switch (p)
        {
            case 0:
            wp[p]=vp[p]*r[2].comp_inv();
                break;
            case 1:
            wp[p]=(vp[p]-wp[p-1]*r[0]-wp[p])*r[2].comp_inv();
                break;
            default:
            wp[p]=(vp[p]-wp[p-1]*r[0]-wp[p]*r[1])*r[2].comp_inv();
                break;
        }
        X+=wp[p]*b;
        p++;
    }
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
    MAX_ITER=300;
}
