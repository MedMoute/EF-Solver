#include "solver.h"
#include "utils.h"

#include <QDir>
#include <QDateTime>

//Sparse Linear Algebra solver system for ComplexVectoriel (C3) Symetrical Matrices
//
//Implemented solvers:
//Iterative Methods
// -Conjugate Gradient
// -MinRES
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
        delete Awp;
        delete wp;
        break;
    case MinRES:
        wp = new VectorC3;
        vp = new VectorC3;
        result = Solve_MinRES();
        delete vp;
        delete wp;
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
        //delete wp,Awp;

        break;
    case MinRES:
        wp = new VectorC3;
        vp = new VectorC3;
        result = Solve_MinRES();
        //delete wp,vp;
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
}



int Solver::Solve_GC()
{
    cout<<"Starting resolution of the problem w/ the Conjugate Gradient Method"<<endl;
    int p=0;
    //GC Initialization
    X=B; //We assume Xo=b
    A->MatVectMul(X,res); //res=AXo
    res-=B; // res= AXo-b
    cout<<"Iter 0 : Err : "<<this->computeError()<<endl;


    w=res;  // wo=res
    //GC iteration
    while (p<MAX_ITER) //Max iteration control
    {
        A->MatVectMul(w,v); //v=Aw
        rho=-(res,w);
        rho/=(v,w);
        X+=w*rho;
        res+=v*rho;
        if((res.Re_norm()/B.Re_norm())<=(eps*eps))
        {
            util::print_separator();
            cout<<"The solver converged."<<endl;
            break;
        }
        p++;
        g=-(res,v)/(v,w);
        w*=g;
        w+=res; //w=res+g*w
        cout<<"Iter "<<p<<" : Err : "<<this->computeError()<<endl;
    }
    //TMP data to grab if the program crashes/cannot recover after quitting the solver
    util::print_separator();
    string filepath;
    filepath = (QDir::currentPath()+"/MINRESsolv_out.tmp").toStdString();
    cout<<"Exporting Mesh File as "<<filepath<<"."<<endl;
    std::ofstream FILE;
    int j;
    FILE.open(filepath.c_str(), std::ios::out);
    if (FILE.fail())
    {
        std::cout<<"Erreur lors de l'ouverture de FILE"<<std::endl;
        return 1;
    }
    else // Working out FILE : dumping output
    {
        FILE<<"#File Generated on "<<QDateTime::currentDateTime().toString().toStdString()<<" by EF-Solver.Solve_GC"<<endl;
    }
    for (j=0;j<B.size();j++)
    {
        FILE<<X.X[0][j]<<"~"<<X.X[1][j]<<"~"<<X.X[2][j]<<endl;
    }
    FILE<<endl;
    for (j=0;j<B.size();j++)
    {
        FILE<<B.X[0][j]<<"~"<<B.X[1][j]<<"~"<<B.X[2][j]<<endl;
    }
    FILE.close();
}

int Solver::Solve_MinRES()
{
    cout<<"Starting resolution of the problem w/ the MINRES Method"<<endl;
    //Since the number of iterations is bunded by MAX_ITER, we can initialize the

    C3 t[MAX_ITER][3];
    C3 r[MAX_ITER][3];

    C3 n_res[MAX_ITER+1];

    C3 c[MAX_ITER];
    C3 s[MAX_ITER];

    C3 a[MAX_ITER];
    C3 b[MAX_ITER];
    wp[0].init(B.size());
    vp[0].init(B.size());
    vp[1].init(B.size());
    //MinRES Initialisation
    for (int p = 1;p<MAX_ITER;p++)
    {
        wp[p].init(B.size());
        vp[p+1].init(B.size());
        if( p==1 )
        {
            //Init
            //We assume X(0)=B
            X=B;

            A->MatVectMul(X,res);
            res-=B;
            vp[p]=res;
            vp[p]/=vp[p].Re_norms();
            wp[p]=vp[p];
            n_res[p]=res.Re_norms();
        }
        A->MatVectMul(vp[p],w);
        if (p>1) {
            t[p-1][2]=t[p][0];
            w=w-vp[p-1]*t[p][0];
        }

        t[p][1]=(w,vp[p]);
        w=w-vp[p]*t[p][1];
        t[p+1][0]=w.Re_norms();


        vp[p+1]=w*(t[p+1][0]).comp_inv();


        //Apply Givens Rotations
        //If p=2 only apply the First rotation
        if (p>2) {
            r[p-2][2]=-t[p-1][2]*s[p-2];
            a[p-1]=t[p-1][2]*c[p-2];

            r[p-1][1]=a[p-1]*c[p-1] - t[p][1]*s[p-1];
            a[p]=a[p-1]*s[p-1] + t[p][1]*c[p-1];
        }
        else
        {
            if (p==2)
            {
                a[p-1]=t[p-1][2];

                r[p-1][1]=a[p-1]*c[p-1] - t[p][1]*s[p-1];
                a[p]=a[p-1]*s[p-1] + t[p][1]*c[p-1];
            }
            else
                if (p==1)
                    a[p]=t[p][1];
        }
        c[p]=C3(C3(C(1,0),C(1,0),C(1,0))/(C3(C(1,0),C(1,0),C(1,0))+(t[p+1][0]/a[p])*(t[p+1][0]/a[p]))).sqrt_();
        s[p]=-c[p]*t[p+1][0]/a[p];
        r[p][0]=c[p]*a[p] -s[p]*t[p+1][0];
        //disp(S(p).C(p))
        b[p]=-c[p]*n_res[p];
        n_res[p+1]=s[p]*n_res[p] ;
        cout<<"Iter num: "<<p<<" - ||res||"<<n_res[p]<<endl;

        //Compute wp

        if (p>2)
        {
            wp[p]=(vp[p]-wp[p-1]*r[p-1][1]-wp[p-2]*r[p-2][2]);
            wp[p]/=r[p][0];
        }
        else
            if (p==2)
            {
                wp[p]=(vp[p]-wp[p-1]*r[p-1][1]);
                wp[p]/=r[p][0];
            }
            else
                if (p==1)
                {
                    wp[p]=vp[p];
                    wp[p]/=r[p][0];
                }
        //compute xp
        X=X+wp[p]*b[p];
        //    cout<<wp[p].X[0]<<endl;
        if(n_res[p].norm()<eps)
            break;
    }
    //TMP data to grab if the program crashes/cannot recover after quitting the solver
    util::print_separator();
    string filepath;
    filepath = (QDir::currentPath()+"/MINRESsolv_out.tmp").toStdString();
    cout<<"Exporting Mesh File as "<<filepath<<"."<<endl;
    std::ofstream FILE;
    int j;
    FILE.open(filepath.c_str(), std::ios::out);
    if (FILE.fail())
    {
        std::cout<<"Erreur lors de l'ouverture de FILE"<<std::endl;
        return 1;
    }
    else // Working out FILE : dumping output
    {
        FILE<<"#File Generated on "<<QDateTime::currentDateTime().toString().toStdString()<<" by EF-Solver.Solve_MinRES"<<endl;
    }
    for (j=0;j<B.size();j++)
    {
        FILE<<X.X[0][j]<<"~"<<X.X[1][j]<<"~"<<X.X[2][j]<<endl;
    }
    FILE<<endl;
    for (j=0;j<B.size();j++)
    {
        FILE<<B.X[0][j]<<"~"<<B.X[1][j]<<"~"<<B.X[2][j]<<endl;
    }
    FILE.close();
}

void Solver::initSolver(int _l)
{
    B.resize(_l);
    X.resize(_l);
    res.resize(_l);
    w.resize(_l);
    v.resize(_l);
    g=C3(0.,0.,0.,0.,0.,0.);
    rho=C3(0.,0.,0.,0.,0.,0.);
    MAX_ITER=100;
}

//DO NOT USE THIS METHOD FOR THE DIRECT SOLVER SINCE THE MATRIX POINTER IS UNSTABLE
//USE THE EXPLICIT METHOD BELOW
R3 Solver::computeError(ErrorNorm norm)
{
    VectorC3 Ax;
    A->MatVectMul(X,Ax);
    VectorC3 ResComp = Ax-B;
    if (norm==ErrorNorm::infinity)
        return(R3(ResComp.X[0].linfty(),ResComp.X[1].linfty(),ResComp.X[2].linfty()));
    else if (norm==ErrorNorm::L1)
        return (R3(ResComp.X[0].l1(),ResComp.X[1].l1(),ResComp.X[2].l1()));
    else if (norm==ErrorNorm::L2)
        return (R3(ResComp.X[0].l2(),ResComp.X[1].l2(),ResComp.X[2].l2()));
    else
    {
        cout<<"Norme de l'erreur mal définie."<<endl;
        return R3(-1,-1,-1);
    }
}

R3 Solver::computeError(MatSparseC3 _A, VectorC3 _B,ErrorNorm norm)
{
    VectorC3 Ax;
    _A.MatVectMul(X,Ax);
    VectorC3 ResComp = Ax-_B;
    if (norm==ErrorNorm::infinity)
        return(R3(ResComp.X[0].linfty(),ResComp.X[1].linfty(),ResComp.X[2].linfty()));
    else if (norm==ErrorNorm::L1)
        return (R3(ResComp.X[0].l1(),ResComp.X[1].l1(),ResComp.X[2].l1()));
    else if (norm==ErrorNorm::L2)
        return (R3(ResComp.X[0].l2(),ResComp.X[1].l2(),ResComp.X[2].l2()));
    else
    {
        cout<<"Norme de l'erreur mal définie."<<endl;
        return R3(-1,-1,-1);
    }
}

void Solver::displaySolution()
{
    cout<<"Solution sur X : "<<endl;
    cout<<X.X[0]<<endl;
    cout<<"Solution sur Y : "<<endl;
    cout<<X.X[1]<<endl;
    cout<<"Solution sur Z : "<<endl;
    cout<<X.X[2]<<endl;
}
