#include "solver.h"
#include "utils.h"

Solver::Solver(FEProblem* _fep, SolverMethod _meth,double _eps):A(_fep->GetA()),eps(_eps)
{
    cout<<"Ajustement des vecteurs a la taille du probleme : "<<_fep->GetRHS().size()<<endl;
    initSolver(_fep->GetRHS().size());
    B=_fep->GetBC()+_fep->GetRHS();
    switch (_meth)
    {
    case GC:
        Solve_GC();
        break;
    case OrthoDir:
        Awp = new VectorC3;
        Solve_Orthodir();
        break;
    case GMRES:
        Solve_GMRES();
    default:
        break;
    }
}

Solver::Solver(MatSparseC3 _mat, VectorC3 _vec, SolverMethod _meth, double _eps):A(&_mat),B(_vec),eps(_eps)
{
    cout<<"Ajustement des vecteurs a la taille du probleme : "<<A->GetSize()<<endl;
    initSolver(A->GetSize());
    switch (_meth)
    {
    case GC:
        Solve_GC();
        break;
    case OrthoDir:
        Awp = new VectorC3;
        Solve_Orthodir();
        break;
    case GMRES:
        Solve_GMRES();
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

    Awp[0].init(B.size());
    Awp[0]=v;
    Awp[0]/=v_n;  //Aw_0=v/||v||
    //    cout<<"Aw1  : "<<Awp[0].X[0]<<endl;
    //    cout<<"res : "<<res.X[0]<<endl;
    rho=-(res,Awp[0]);

    res+=Awp[0]*rho;
    X+=w*rho;
    cout<<X.X[0]<<endl;


    j=0;
    while (res.Re_norm()>eps && j<50)
    {
        w=Awp[j]; //w=Aw_p-1
        A->MatVectMul(w,v); //v=Aw
        //cout<<"v ["<<j<<"] : "<<v.X[0]<<endl;
            //Debugged up to here

        //Orhogonalisation de la direction de descente par rapport aux directions précedentes
        for(p=0;p<=j;p++)
        {
            C3 alpha=(v,Awp[p]); //a_p = (v.Aw_p)
            cout<<"Alpha : "<<alpha.X_()<<endl;
            w-=v*alpha; // w = w-a_p*v_p
            v-=Awp[p]*alpha; // v = v-a_p*Aw_p
            cout<<"v ["<<j<<"] : "<<v.X[0]<<endl;
            cout<<"w ["<<j<<"] : "<<w.X[0]<<endl;
        }
        //cout<<"w  : "<<w.X[0]<<endl;
        //Orthonormalisation
        R3 v_n=v.Re_norms();
        w/=v_n;  // w_p+1= w /||v||

        Awp[j+1].init(B.size());
        Awp[j+1]=v;
        Awp[j+1]/=v_n; // Aw_p+1= v /||v||
        //cout<<"v_n = "<<v_n.X_()<<"-"<<v_n.Y_()<<"-"<<v_n.Z_()<<endl;

        rho=-(res,Awp[j+1]);
        // Debug : Verif de l'AtA- orthogonalisation des résidus
        //                for(p=0;p<=j+1;p++)
        //                {
        //                    C3 test= (Awp[p],Awp[j+1]);
        //                    cout<<"Ortho des directions de la base a l'etape "<<j+1<<" vs "<<p<<" : "<<test.norm()<<endl;
        //                }
        //        cout<<rho.X_()<< " - "<<rho.Y_()<< " - "<<rho.Z_()<<endl;
        X+=w*rho;
        //cout<<X.X[0]<<endl;

        res+=Awp[j+1]*rho;
        /*
         *Debug (check if the Krylov basis vectors are A^tA-orthogonaux)
         *for(p=1;p<=j;p++)
        {
            C3 test= (Awp[p],res);
            cout<<"Ortho du residu a l'etape "<<j+1<<" vs "<<p<<" : "<<test.norm()<<endl;
        }
        */

        cout<<"Iter "<<j+1<<": ||res|| ="<<res.Re_norm()<<endl;
        j++;
    }
    cout<<"Solution sur X : "<<endl;
    cout<<X.X[0]<<endl;
}

int Solver::Solve_GC()
{
    cout<<"Starting resolution of the problem w/ the Conjugate Gradient Method"<<endl;
    int p=0;
    //GC Initialieation
    X=B; //We assume Xo=b
    A->MatVectMul(X,res); //res=AXo
    res-=B; // res= AXo-b
    cout<<res.X[0]<<endl;
    cout<<"Iter 0 : ||res|| = "<<res.Re_norms().X_()<<endl;


    w=res;  // wo=res
    //GC iteration
    while (p<50) //Max iteration control
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
            cout<<X.X[2]<<endl;            break;
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
}
