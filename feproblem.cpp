#include "feproblem.h"
#include "C3.h"
#include "utils.h"
#include "ui_parametersdialog.h"

#include <chrono>

FEProblem::FEProblem(Maillage3D* _maillage,RHS* _rhs,BC* _bc,int _op,bool _unit):maillage(_maillage),op(_op),unitary(_unit)
{
    A = new MatSparseC3;
    util::print_separator();
    A->EditProfileFromMesh(*maillage);

    auto start_time = std::chrono::system_clock::now();
    //
    //Filling of vector G from BC data
    //We use the gauss integration formula for quadratic functions over a simplex
    //since the node function v_i(X) is piecewise linear, we'll have exact integration
    //value for piecewise affines boundary conditions,which is precisely the case
    //for imported data from the mesh geometry. If we have an analytic expression
    // we'll use it instead of the interpolated data on the mesh.

    cout<<"Filling of vector G from BC data"<<endl;

    map<int,R3>* b_n = maillage->GetBoundaryNodesMap();
    map<int,R3>* n = maillage->GetNodesMap();
    map<int,R3>::iterator it_n;
    map<int,R3>::iterator it_b_n;
    map<int,set<int>*>::iterator it_map;
    set<int>* set_b_n;
    set<int>::iterator it_set;
    set<int>* set_n;

    int max_n_idx =n->rbegin()->first;
    int i;
    for (i=0;i<3;i++)
    {
        F.X[i].resize(max_n_idx);
        G.X[i].resize(max_n_idx);
    }
    if(_bc!=0)
    {
        int i;
        int j=0;
        for(it_b_n=b_n->begin();it_b_n!=b_n->end();++it_b_n)
        {//Loop on the boundary nodes
            int idx=it_b_n->first;
            //Fetch the faces which are adjacent to that node
            it_map=maillage->GetNeighboursFacesMap()->find(idx);
            if(it_map !=maillage->GetNeighboursFacesMap()->end())
            {
                set_b_n=it_map->second;
                for(it_set=set_b_n->begin();it_set!=set_b_n->end();++it_set)
                {
                    C3 res =CalcIntOnTriangle(maillage,_bc,idx,*it_set);
                    for (i=0;i<3;i++)
                    {
                        (G.X[i])[idx-1]+=C(res[i]);
                    }
                }
            }
            j++;
            if (j%100==0)
            {
                cout<<"["<<int(j)*100/b_n->size()<<"%]"<<"Computing Gauss Integration of triangle #"<<j<<"/"<<b_n->size()<<" ."<<endl;
            }
        }
        cout <<"Boundary conditions fully determined.Proceeding to next step."<<endl;
    }
    util::print_separator();

    cout<<"Filling of vector F from RHS data"<<endl;
    //Filling of vector F from RHS data

    if(_rhs!=0)
    {
        int i;
        int j=0;

        for(it_n=n->begin();it_n!=n->end();++it_n)
        {//Loop on the boundary nodes
            int idx=it_n->first;
            //Fetch the elements which are adjacent to that node
            it_map=maillage->GetNeighboursElementsMap()->find(idx);
            if(it_map !=maillage->GetNeighboursElementsMap()->end())
            {
                set_n=it_map->second;
                for(it_set=set_n->begin();it_set!=set_n->end();++it_set)
                {
                    C3 res = CalcIntOnTetrahedron(maillage,_rhs,idx,*it_set);
                    for (i=0;i<3;i++)
                    {
                        (F.X[i])[idx-1]+=C(res[i]);
                    }
                }
            }
            j++;
            if (j%100==0)
            {
                cout<<"["<<int(i)*100/n->size()<<"%]"<<"Computing Gauss Integration of tetrahedron #"<<i<<"/"<<n->size()<<"."<<endl;
            }
        }
        cout <<"RHS Member fully determined.Proceeding to next step."<<endl;

    }
    auto end_time = std::chrono::system_clock::now();

    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Calcul du second membre et des conditions aux bord matricielles : effectue en "<<elapsed_seconds<<"s."<<std::endl;
    int p=1;
    util::print_separator();
    cout<<"Remplissage de la matrice sparse de l'opérateur"<<endl;
    //Remarque : La matrice A de l'opérateur n'est jamais construite densément
    //on construit d'abord son profil creux, puis on rempit les données par recherche d'indices dans le profil creux
    switch (op)
    {
    case 0 :
    {
        cout<<"Construction de l'operateur pour l'equiation de Poisson"<<endl;
        if(unitary)
            BuildPoissonOperator();
        else//Parametrized problem
        {
            //Open dialog.
            dialog=new ParametersDialog(maillage->GetPartitionsMap());
            dialog->show();
            QObject::connect(dialog->GetUI()->pushButton_ok,SIGNAL(clicked()),this,SLOT(ProcessParamData()));
        }

        break;
    }

    case 1 :
    {
        cout<<"Construction de l'operateur pour l'equation de Helmholtz"<<endl;
        if(unitary)
            BuildHelmholtzOperator();
        else//Parametrized problem
        {
            //Open dialog.
            dialog=new ParametersDialog(maillage->GetPartitionsMap());
            dialog->show();
            QObject::connect(dialog->GetUI()->pushButton_ok,SIGNAL(clicked()),this,SLOT(ProcessParamData()));
        }

        break;
    }
    case 2 :
    {
        cout<<"Construction de l'operateur pour l'equiation de la chaleur stationnaire"<<endl;
        if(unitary)
            BuildHeatOperator();
        else//Parametrized problem
        {
            //Open dialog.
            dialog=new ParametersDialog(maillage->GetPartitionsMap());
            dialog->show();
            QObject::connect(dialog->GetUI()->pushButton_ok,SIGNAL(clicked()),this,SLOT(ProcessParamData()));
        }
        break;
    }
    default:
        break;
    }

}

C3 FEProblem::CalcIntOnTriangle(Maillage3D* _maillage,BC* _bc,int _idx,int _tri_idx)
{
    map<int,tuple<int,int,int>>::iterator face =_maillage->GetBoundaryFacesMap()->find(_tri_idx);
    map<int,R3>* nodes = _maillage->GetNodesMap();
    int tri_n_idx[3];
    int local_idx=-1;
    R pos_x_n[3];
    R pos_y_n[3];
    R pos_z_n[3];

    R pos_x_g[3];
    R pos_y_g[3];
    R pos_z_g[3];

    R val_v_g[3];

    C tmp_val_g[9];

    C3 val_g[3];

    C3 res;
    int i,j;
    if(face!=_maillage->GetBoundaryFacesMap()->end())
    {
        tri_n_idx[0]=get<0>(face->second);
        tri_n_idx[1]=get<1>(face->second);
        tri_n_idx[2]=get<2>(face->second);
        for (i=0;i<3;i++)
        {
            R3 pos_n = nodes->find(tri_n_idx[i])->second;
            pos_x_n[i]=pos_n.X_();
            pos_y_n[i]=pos_n.Y_();
            pos_z_n[i]=pos_n.Z_();
            if(tri_n_idx[i]==_idx)
            {
                local_idx=i;
            }
        }
        //On effectue l'integration sur les points de gauss suivants:
        //(1/6,1/6) - (2/3,1/6) - (1/6,2/3)
        for (i=0;i<3;i++)
        {
            pos_x_g[i]=(1./6)*pos_x_n[i]+(1./6)*pos_x_n[(i+1)%3]+(2./3)*pos_x_n[(i+2)%3];
            pos_y_g[i]=(1./6)*pos_y_n[i]+(1./6)*pos_y_n[(i+1)%3]+(2./3)*pos_y_n[(i+2)%3];
            pos_z_g[i]=(1./6)*pos_z_n[i]+(1./6)*pos_z_n[(i+1)%3]+(2./3)*pos_z_n[(i+2)%3];
        }
        string* bc_expr = _bc->GetExpr();
        //Calcul des valeurs en pos_g[i] à partir de l'expression de bc si possible

        if(bc_expr[0].size()>0 && bc_expr[1].size()>0 && bc_expr[2].size()>0)
        {
            unsigned int i;
            ParserX  p(pckALL_COMPLEX);
            Value xVal;
            Value yVal;
            Value zVal;
            p.DefineVar("x",  Variable(&xVal));
            p.DefineVar("y",  Variable(&yVal));
            p.DefineVar("z",  Variable(&zVal));


            for (i=0;i<3;i++)
            {
                // cout<<pos_x_g[i]<<pos_y_g[i]<<pos_z_g[i]<<endl;

                xVal=pos_x_g[i];
                yVal=pos_y_g[i];
                zVal=pos_z_g[i];

                for (j=0;j<3;j++)
                {
                    p.SetExpr(bc_expr[j]);
                    tmp_val_g[3*i+j] =p.Eval().GetComplex();
                    // cout<<bc_expr[j]<<" @ "<<xVal<<"/"<<yVal<<"/"<<zVal<<": "<<p.Eval().GetComplex()<<endl;
                }
                val_g[i]=C3(tmp_val_g[3*i],tmp_val_g[3*i+1],tmp_val_g[3*i+2]);
            }
        }
        //sinon interpolation lineaire sur l'array de données
        else
        {
            vector<C3>* data =_bc->GetBC();
            C3 val_n[3];
            for(i=0;i<3;i++)
            {
                val_n[i]=data->at(tri_n_idx[i]);
            }
            for(i=0;i<3;i++)
            {
                val_g[i]=val_n[i]*(1./6)+val_n[(i+1)%3]*(1./6)+val_n[(i+2)%3]*(2./3);
            }
        }
        for (i=0;i<3;i++)
        {
            if (local_idx==i)
            {
                val_v_g[i]=2./3;
                val_v_g[(i+1)%3]=1./6;
                val_v_g[(i+2)%3]=1./6;
            }
        }
        res =(val_g[0]*val_v_g[0]+val_g[1]*val_v_g[1]+val_g[2]*val_v_g[2])*(1./6);
        //cout<<"Node "<<_idx<<"+="<<res.comp(1)<<";"<<res.comp(2)<<";"<<res.comp(3)<<";"<<endl;
    }
    return res;
}

C3 FEProblem::CalcIntOnTetrahedron(Maillage3D* _maillage,RHS* _rhs,int _idx,int _tetra_idx)
{
    map<int,tuple<int,int,int,int>>::iterator elem =_maillage->GetElementsMap()->find(_tetra_idx);
    map<int,R3>* nodes = _maillage->GetNodesMap();
    int tri_n_idx[4];
    int local_idx=-1;
    R pos_x_n[4];
    R pos_y_n[4];
    R pos_z_n[4];

    R pos_x_g[4];
    R pos_y_g[4];
    R pos_z_g[4];

    R val_v_g[4];

    C tmp_val_g[12];

    C3 val_g[3];

    C3 res;
    int i,j;
    if(elem!=_maillage->GetElementsMap()->end())
    {
        tri_n_idx[0]=get<0>(elem->second);
        tri_n_idx[1]=get<1>(elem->second);
        tri_n_idx[2]=get<2>(elem->second);
        tri_n_idx[4]=get<3>(elem->second);

        for (i=0;i<4;i++)
        {
            R3 pos_n = nodes->find(tri_n_idx[i])->second;
            pos_x_n[i]=pos_n.X_();
            pos_y_n[i]=pos_n.Y_();
            pos_z_n[i]=pos_n.Z_();
            if(tri_n_idx[i]==_idx)
            {
                local_idx=i;
            }
        }
        //On effectue l'integration sur les points de gauss suivants:
        //w1 = 0.25 ;g1(0.1381966011250105,0.1381966011250105,0.1381966011250105,0.5854101966249685)
        //w2 = 0.25 ;g2(0.1381966011250105,0.1381966011250105,0.5854101966249685,0.1381966011250105)
        //w3 = 0.25 ;g3(0.1381966011250105,0.5854101966249685,0.1381966011250105,0.1381966011250105)
        //w4 = 0.25 ;g4(,0.58541019662496850.1381966011250105,0.1381966011250105,0.1381966011250105)
        for (i=0;i<4;i++)
        {
            pos_x_g[i]=(0.1381966011250105)*pos_x_n[i]+(0.1381966011250105)*pos_x_n[(i+1)%4]+(0.1381966011250105)*pos_x_n[(i+2)%4]+(0.5854101966249685)*pos_x_n[(i+3)%4];
            pos_y_g[i]=(0.1381966011250105)*pos_y_n[i]+(0.1381966011250105)*pos_y_n[(i+1)%4]+(0.1381966011250105)*pos_y_n[(i+2)%4]+(0.5854101966249685)*pos_x_n[(i+3)%4];;
            pos_z_g[i]=(0.1381966011250105)*pos_z_n[i]+(0.1381966011250105)*pos_z_n[(i+1)%4]+(0.1381966011250105)*pos_z_n[(i+2)%4]+(0.5854101966249685)*pos_x_n[(i+3)%4];;
        }
        string* rhs_expr = _rhs->GetExpr();
        //Calcul des valeurs en pos_g[i] à partir de l'expression de bc si possible

        if(rhs_expr[0].size()>0 && rhs_expr[1].size()>0 && rhs_expr[2].size()>0)
        {
            unsigned int i;
            ParserX  p(pckALL_COMPLEX);
            Value xVal;
            Value yVal;
            Value zVal;
            p.DefineVar("x",  Variable(&xVal));
            p.DefineVar("y",  Variable(&yVal));
            p.DefineVar("z",  Variable(&zVal));


            for (i=0;i<4;i++)
            {
                // cout<<pos_x_g[i]<<pos_y_g[i]<<pos_z_g[i]<<endl;

                xVal=pos_x_g[i];
                yVal=pos_y_g[i];
                zVal=pos_z_g[i];

                for (j=0;j<3;j++)
                {
                    p.SetExpr(rhs_expr[j]);
                    tmp_val_g[3*i+j] =p.Eval().GetComplex();
                    // cout<<bc_expr[j]<<" @ "<<xVal<<"/"<<yVal<<"/"<<zVal<<": "<<p.Eval().GetComplex()<<endl;
                }
                val_g[i]=C3(tmp_val_g[3*i],tmp_val_g[3*i+1],tmp_val_g[3*i+2]);
            }
        }
        //sinon interpolation lineaire sur l'array de données
        else
        {
            vector<C3>* data =_rhs->GetRHS();
            C3 val_n[3];
            for(i=0;i<4;i++)
            {
                val_n[i]=data->at(tri_n_idx[i]);
            }
            for(i=0;i<4;i++)
            {
                val_g[i]=val_n[i]*(0.1381966011250105)+val_n[(i+1)%4]*(0.1381966011250105)+val_n[(i+2)%4]*(0.1381966011250105)+val_n[(i+3)%4]*(0.5854101966249685);
            }
        }
        for (i=0;i<3;i++)
        {
            if (local_idx==i)
            {
                val_v_g[i]=0.5854101966249685;
                val_v_g[(i+1)%4]=0.1381966011250105;
                val_v_g[(i+2)%4]=0.1381966011250105;
                val_v_g[(i+3)%4]=0.1381966011250105;
            }
        }
        res =(val_g[0]*val_v_g[0]+val_g[1]*val_v_g[1]+val_g[2]*val_v_g[2]+val_g[3]*val_v_g[3])*(1./(6.*4.));
        //cout<<"Node "<<_idx<<"+="<<res.comp(1)<<";"<<res.comp(2)<<";"<<res.comp(3)<<";"<<endl;
    }
    return res;
}


void FEProblem::DisplayElementaryMassMatrix(int n_elem)
{
    util::print_separator();
    cout<<"Affichage de la matrice de masse de l'element "<<n_elem<<endl<<endl;
    map<int,tuple<int,int,int,int>>* m_ptr = maillage->GetElementsMap();
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    m_it = m_ptr->find(n_elem);
    int nodes[4];
    R3 pos[4];
    C Mat[16];
    int i,j;
    if(m_it!=m_ptr->end())//Check is the element exists
    {
        tuple<int,int,int,int> elem = m_it->second;

        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            cout<<nodes[i]<<" : ";
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
                cout<<pos[i].X_()<<"-"<<pos[i].Y_()<<"-"<<pos[i].Z_()<<endl;
            }
            else
            {
                cout<<"Il n'y a pas de noeud correspondant a la valeur recherchee."<<endl;
            }
        }
    }
    else
    {
        cout<<"Il n'y a pas d'element correspondant a la valeur recherchee."<<endl;
    }
    double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
    cout<<"Mesure signee de l'element: "<<V<<endl;
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            Mat[4*i+j]=abs(V)*(1+util::kron(i,j))/20;
            cout<<Mat[4*i+j]<<" ";
        }
        cout<<endl;
    }
}

void FEProblem::DisplayElementaryStiffnessMatrix(int n_elem)
{
    util::print_separator();
    cout<<"Affichage de la matrice de rigidite de l'element "<<n_elem<<endl<<endl;
    map<int,tuple<int,int,int,int>>* m_ptr = maillage->GetElementsMap();
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    m_it = m_ptr->find(n_elem);
    int nodes[4];
    R3 pos[4];
    double a[4];
    double b[4];
    double c[4];
    C Mat[16];
    int i,j;
    if(m_it!=m_ptr->end())//Check is the element exists
    {
        tuple<int,int,int,int> elem = m_it->second;

        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            cout<<nodes[i]<<" : ";
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check is the node exists
            {
                pos[i]=m_pos_it->second;
                cout<<pos[i].X_()<<"-"<<pos[i].Y_()<<"-"<<pos[i].Z_()<<endl;
            }
            else
            {
                cout<<"Il n'y a pas de noeud correspondant a la valeur recherchee."<<endl;
            }
        }
    }
    else
    {
        cout<<"Il n'y a pas d'element correspondant a la valeur recherchee."<<endl;
    }
    double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
    cout<<"Mesure signee de l'element: "<<V<<endl;
    for(i=0;i<4;i++)
    {
        a[i]=(pos[(i+1)%4].Y_()-pos[(i+3)%4].Y_())*pos[(i+2)%4].Z_();
        b[i]=(pos[(i+3)%4].X_()*pos[(i+1)%4].Z_())-(pos[(i+1)%4].X_()*pos[(i+3)%4].Z_());
        c[i]=(pos[(i+3)%4].Y_()-pos[(i+1)%4].Y_())*pos[(i+2)%4].X_();
    }
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            Mat[4*i+j]=(a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.);
            cout<<Mat[4*i+j]<<" ";
        }
        cout<<endl;
    }
}
//Construction de l'opérateur de Poisson en
//utilisant la matrice de masse
void FEProblem::BuildPoissonOperator()
{
    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];
    R3 pos[4];
    int i,j;
    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }

        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C(abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C(abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C(abs(V)*(1+util::kron(i,j))/20,0);
            }
        }
        //Debug :
        //DisplayElementaryMassMatrix(n_elem);
    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de Laplace : effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::BuildParameterPoissonOperator(double lambda)
{
    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];
    R3 pos[4];
    int i,j;
    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }

        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C(lambda*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C(lambda*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C(lambda*abs(V)*(1+util::kron(i,j))/20,0);
            }
        }
        //Debug :
        //DisplayElementaryMassMatrix(n_elem);
    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de Laplace : effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::BuildPartitionnedParameterPoissonOperator(std::map<int,double> lambdas,std::map<int,int>* partition_data)
{
    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];
    R3 pos[4];
    int i,j;
    int part;
    double elem_data;
    map<int,int>::iterator part_it;
    map<int,double>::iterator part_val_it;

    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }

        //Get the partition index of the current element :
        part_it= partition_data->find(n_elem);
        part=(*part_it).second;
        //find the data corresponding to the partition

        for(part_val_it=lambdas.begin();part_val_it!=lambdas.end();++part_val_it)
        {
            if((*part_val_it).first ==part)
                elem_data=(*part_val_it).second;
        }

        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C(elem_data*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C(elem_data*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C(elem_data*abs(V)*(1+util::kron(i,j))/20,0);
            }
        }
        //Debug :
        //DisplayElementaryMassMatrix(n_elem);
    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de Laplace : effectue en "<<elapsed_seconds<<"s."<<std::endl;
}
//Construction de l'opérateur de la Chaleur Stationnaire en
//utilisant la matrice de rigidité du système
void FEProblem::BuildHeatOperator()
{

    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];

    double a[4];
    double b[4];
    double c[4];
    R3 pos[4];
    int i,j;
    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }

        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            a[i]=(pos[(i+1)%4].Y_()-pos[(i+3)%4].Y_())*pos[(i+2)%4].Z_();
            b[i]=(pos[(i+3)%4].X_()*pos[(i+1)%4].Z_())-(pos[(i+1)%4].X_()*pos[(i+3)%4].Z_());
            c[i]=(pos[(i+3)%4].Y_()-pos[(i+1)%4].Y_())*pos[(i+2)%4].X_();
        }
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
            }
        }
    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de la Chaleur Stationnaire : effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::BuildParameterHeatOperator(double lambda)
{

    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];

    double a[4];
    double b[4];
    double c[4];
    R3 pos[4];
    int i,j;
    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }


        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            a[i]=(pos[(i+1)%4].Y_()-pos[(i+3)%4].Y_())*pos[(i+2)%4].Z_();
            b[i]=(pos[(i+3)%4].X_()*pos[(i+1)%4].Z_())-(pos[(i+1)%4].X_()*pos[(i+3)%4].Z_());
            c[i]=(pos[(i+3)%4].Y_()-pos[(i+1)%4].Y_())*pos[(i+2)%4].X_();
        }
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C(lambda*(a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C(lambda*(a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C(lambda*(a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
            }
        }
    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de la Chaleur Stationnaire : effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::BuildPartitionnedParameterHeatOperator(std::map<int, double> lambdas,std::map<int,int>* partition_data)
{
    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];

    double a[4];
    double b[4];
    double c[4];
    R3 pos[4];
    int i,j;
    int part;
    double elem_data;
    map<int,int>::iterator part_it;
    map<int,double>::iterator part_val_it;

    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }
        //Get the partition index of the current element :
        part_it= partition_data->find(n_elem);
        part=(*part_it).second;
        //find the data corresponding to the partition

        for(part_val_it=lambdas.begin();part_val_it!=lambdas.end();++part_val_it)
        {
            if((*part_val_it).first ==part)
                elem_data=(*part_val_it).second;
        }
        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            a[i]=(pos[(i+1)%4].Y_()-pos[(i+3)%4].Y_())*pos[(i+2)%4].Z_();
            b[i]=(pos[(i+3)%4].X_()*pos[(i+1)%4].Z_())-(pos[(i+1)%4].X_()*pos[(i+3)%4].Z_());
            c[i]=(pos[(i+3)%4].Y_()-pos[(i+1)%4].Y_())*pos[(i+2)%4].X_();
        }
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C(elem_data*(a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C(elem_data*(a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C(elem_data*(a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
            }
        }
    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de la Chaleur Stationnaire : effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::BuildHelmholtzOperator()
{
    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];
    double a[4];
    double b[4];
    double c[4];
    R3 pos[4];
    int i,j;
    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }

        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C(abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C(abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C(abs(V)*(1+util::kron(i,j))/20,0);
            }
        }
        //Add the stiffness part
        for(i=0;i<4;i++)
        {
            a[i]=(pos[(i+1)%4].Y_()-pos[(i+3)%4].Y_())*pos[(i+2)%4].Z_();
            b[i]=(pos[(i+3)%4].X_()*pos[(i+1)%4].Z_())-(pos[(i+1)%4].X_()*pos[(i+3)%4].Z_());
            c[i]=(pos[(i+3)%4].Y_()-pos[(i+1)%4].Y_())*pos[(i+2)%4].X_();
        }
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
            }
        }

    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de Helmholz: effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::BuildParameterHelmholtzOperator(double omega)
{
    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];
    double a[4];
    double b[4];
    double c[4];
    R3 pos[4];
    int i,j;
    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }

        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C(omega*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C(omega*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C(omega*abs(V)*(1+util::kron(i,j))/20,0);
            }
        }
        //Add the stiffness part
        for(i=0;i<4;i++)
        {
            a[i]=(pos[(i+1)%4].Y_()-pos[(i+3)%4].Y_())*pos[(i+2)%4].Z_();
            b[i]=(pos[(i+3)%4].X_()*pos[(i+1)%4].Z_())-(pos[(i+1)%4].X_()*pos[(i+3)%4].Z_());
            c[i]=(pos[(i+3)%4].Y_()-pos[(i+1)%4].Y_())*pos[(i+2)%4].X_();
        }
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
            }
        }

    }
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de Helmholz: effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::BuildPartitionnedParameterHelmholtzOperator(std::map<int,double> omegas,std::map<int,int>* partition_data)
{
    //Timer start
    auto start_time = std::chrono::system_clock::now();

    map<int,tuple<int,int,int,int>>::iterator m_it;
    map<int,R3>::iterator m_pos_it;
    map<int,R3>* m_pos_ptr = maillage->GetNodesMap();
    int n_elem;
    int nodes[4];
    int nodes_tri[3];
    double a[4];
    double b[4];
    double c[4];
    R3 pos_tri[3];
    R3 pos[4];
    int i,j;

    int part;
    double elem_data;
    map<int,int>::iterator part_it;
    map<int,double>::iterator part_val_it;

    for(m_it=maillage->GetElementsMap()->begin();m_it!=maillage->GetElementsMap()->end();++m_it)
    {
        n_elem=m_it->first;
        tuple<int,int,int,int> elem = m_it->second;
        nodes[0]=get<0>(elem);
        nodes[1]=get<1>(elem);
        nodes[2]=get<2>(elem);
        nodes[3]=get<3>(elem);

        for  (i=0;i<4;i++)
        {
            m_pos_it = m_pos_ptr->find(nodes[i]);
            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
            {
                pos[i]=m_pos_it->second;
            }
        }

        //Get the partition index of the current element :
        part_it= partition_data->find(n_elem);
        part=(*part_it).second;
        //find the data corresponding to the partition

        for(part_val_it=omegas.begin();part_val_it!=omegas.end();++part_val_it)
        {
            if((*part_val_it).first ==part)
                elem_data=(*part_val_it).second;
        }

        double V=util::simplex3Mesure(pos[0],pos[1],pos[2],pos[3]);
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]-=C(elem_data*elem_data*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]-=C(elem_data*elem_data*abs(V)*(1+util::kron(i,j))/20,0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]-=C(elem_data*elem_data*abs(V)*(1+util::kron(i,j))/20,0);
            }
        }
        //Add the stiffness part
        for(i=0;i<4;i++)
        {
            a[i]=(pos[(i+1)%4].Y_()-pos[(i+3)%4].Y_())*pos[(i+2)%4].Z_();
            b[i]=(pos[(i+3)%4].X_()*pos[(i+1)%4].Z_())-(pos[(i+1)%4].X_()*pos[(i+3)%4].Z_());
            c[i]=(pos[(i+3)%4].Y_()-pos[(i+1)%4].Y_())*pos[(i+2)%4].X_();
        }
        for(i=0;i<4;i++)
        {
            for(j=0;j<4;j++)
            {
                A->GetData().X[0][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[1][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
                A->GetData().X[2][A->IdxSearch(nodes[i],nodes[j])]+=C((a[i]*a[j]+b[i]*b[j]+c[i]*c[j])/(abs(V)*36.),0);
            }
        }
    }
    //Add the interface part: we use another loop on the smaller map : triangles_interfaces et triangles_interfaces_nodes
    map<int,pair<int,int>>::iterator it;
    int tri;
    int part_i,part_j;
    std::map<int,double>::iterator m_part_it ;
    std::map<int,double>::iterator m_o_part_it;
    std::map<int,pair<int,int>>::iterator t_it;
    std::map<int,tuple<int,int,int>>::iterator t_n_it;
    for (i=0;i<maillage->GetPartitionsMap()->size();i++)
    {
        part_i=maillage->GetPartitionsMap()->at(i);
        for (j=0;j<maillage->GetPartitionsMap()->size();j++)
        {
            part_j=maillage->GetPartitionsMap()->at(j);
            if (part_i!=part_j)
            {
                // Interface part_i/part_j

                m_part_it = omegas.find(part_i);
                m_o_part_it = omegas.find(part_j);
                elem_data=(*m_part_it).second-(*m_o_part_it).second;

                //Parcours des éléments triangles sur l'interface en cours
                for (t_it = maillage->GetTrianglesInterfaceMap()->begin(); t_it != maillage->GetTrianglesInterfaceMap()->end(); ++t_it )
                {
                    if (t_it->second == pair<int,int>(part_i,part_j))
                    {
                        tri = t_it->first;
                        //We now have obtained a triangle that is on the interface between part_i and part_j.
                        // We'll get the positions of its nodes

                        t_n_it = maillage->GetTrianglesNodesInterfaceMap()->find(tri);
                        tuple<int,int,int> tri = t_n_it->second;
                        nodes_tri[0]=get<0>(tri);
                        nodes_tri[1]=get<1>(tri);
                        nodes_tri[2]=get<2>(tri);

                        for  (i=0;i<3;i++)
                        {
                            m_pos_it = m_pos_ptr->find(nodes_tri[i]);
                            if(m_pos_it!=m_pos_ptr->end())//Check if the node exists
                            {
                                pos_tri[i]=m_pos_it->second;
                            }
                        }
                        double a=util::simplex2Mesure(pos[0],pos[1],pos[2]);
                        for(i=0;i<3;i++)
                        {
                            for(j=0;j<3;j++)
                            {
                                A->GetData().X[0][A->IdxSearch(nodes_tri[i],nodes_tri[j])]-=C(0,elem_data*abs(a)/8);
                                A->GetData().X[1][A->IdxSearch(nodes_tri[i],nodes_tri[j])]-=C(0,elem_data*abs(a)/8);
                                A->GetData().X[2][A->IdxSearch(nodes_tri[i],nodes_tri[j])]-=C(0,elem_data*abs(a)/8);
                            }
                        }
                    }
                }
            }
        }
        cout<<endl;
    }
    util::print_separator();
    auto end_time = std::chrono::system_clock::now();
    auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
    std::cout<<"Constuction de l'operateur de Helmholz: effectue en "<<elapsed_seconds<<"s."<<std::endl;
}

void FEProblem::ProcessParamData()
{
    //
    bool global=dialog->GetUI()->radioButton_global->isChecked();
    switch (op)
    {
    case 0 :
    {
        if (global)
            BuildParameterPoissonOperator(dialog->GetUI()->doubleSpinBox->value());
        else
            BuildPartitionnedParameterPoissonOperator(dialog->GetData(),maillage->GetPartitionDataMap());

    }

    case 1 :
    {
        if (global)
            BuildParameterHelmholtzOperator(dialog->GetUI()->doubleSpinBox->value());
        else
            BuildPartitionnedParameterHelmholtzOperator(dialog->GetData(),maillage->GetPartitionDataMap());
    }
    case 2 :
    {
        if (global)
            BuildParameterHeatOperator(dialog->GetUI()->doubleSpinBox->value());
        else
            BuildPartitionnedParameterHeatOperator(dialog->GetData(),maillage->GetPartitionDataMap());
    }
    default:
        break;
    }

}
