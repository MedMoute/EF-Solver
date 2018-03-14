#include "MatSparseC3.h"

VectorC3::VectorC3(ifstream &fd)
{
    //loading a sparse Matrix from a txt File
    // no consistency checks are made for now
    string line;
    int i,data_size;
    if(fd.is_open() ) {
        // on lit les lignes une par une et on avise en fonction des mots cle
        while ( fd.good() ) {
            getline (fd,line);
            // on regarde si on a un mot cle
            if( line[0] == '$' ) {

                if( line.compare("$VData")==0 )
                {
                    //cout<<"VData Processing"<<endl;
                    fd>>data_size;
                    this->resize(data_size);
                    getline(fd,line);
                    for( i=0; i<data_size ; i++)
                    {
                        // on ne sais pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
                        getline(fd,line) ;
                        // on compte les espaces pour savoir combien on a de parametres
                        int nspace=0;
                        for(unsigned int j=0; j!=line.size(); ++j)
                            nspace+=( line.at(j)==' ');
                        int ntmp = nspace+1 ;
                      //  cout<<line<<": Nb espaces = "<<ntmp<<endl;
                        if (ntmp==6)
                        {
                            double* tmp = new double[ntmp] ;
                            // on met les entiers dans les cases
                            size_t spos ;
                            size_t isize ;
                            size_t sdeb = 0 ;
                            int ind = 0 ;
                            do {
                                spos = line.find(' ', sdeb) ;
                                isize = spos - sdeb ;
                                istringstream(line.substr(sdeb, isize)) >> tmp[ind] ;
                                ind ++;
                                sdeb = spos+1 ;
                            } while (spos!=string::npos ) ;
                            // on a toutes les composantes dans tmp
                            X[0][i]=C(tmp[0],tmp[3]);
                            X[1][i]=C(tmp[1],tmp[4]);
                            X[2][i]=C(tmp[2],tmp[5]);
                        //    cout<<X[0][i]<<"  -   "<<X[1][i]<<"   -   "<<X[2][i]<<endl;
                        }
                    }
                }

            }
        }
    }
}

MatSparseC3::MatSparseC3()
{}

MatSparseC3::MatSparseC3(ifstream &fd)
{
    //loading a sparse Matrix from a txt File
    // no consistency checks are made for now
    string line;
    int _size,i,data_size;
    if(fd.is_open() ) {
        // on lit les lignes une par une et on avise en fonction des mots cle
        while ( fd.good() ) {
            getline (fd,line);
            // on regarde si on a un mot cle
            if( line[0] == '$' ) {

                if( line.compare("$PData")==0 )
                {
                    fd>>_size;
                    this->size=_size;
                    P.resize(size);
                    getline(fd,line);
                    for( i=0; i<size ; i++)
                    {
                        // on ne sais pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
                        getline(fd,line) ;
                        // on compte les espaces pour savoir combien on a de parametres
                        int nspace=0;
                        for(unsigned int j=0; j!=line.size(); ++j)
                            nspace+=( line.at(j)==' ');
                        int ntmp = nspace+1 ;
                        if (ntmp=1)
                        {
                            double* tmp = new double[ntmp] ;
                            // on met les entiers dans les cases
                            size_t spos ;
                            size_t isize ;
                            size_t sdeb = 0 ;
                            int ind = 0 ;
                            do {
                                spos = line.find(' ', sdeb) ;
                                isize = spos - sdeb ;
                                istringstream(line.substr(sdeb, isize)) >> tmp[ind] ;
                                ind ++;
                                sdeb = spos+1 ;
                            } while (spos!=string::npos ) ;
                            // on a toutes les composantes dans tmp
                            P[i]=tmp[0];
                        }
                    }
                }

                if( line.compare("$JData")==0 )
                {
                    fd>>data_size;
                    J.resize(data_size);
                    getline(fd,line);
                    for( i=0; i<data_size ; i++)
                    {
                        // on ne sais pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
                        getline(fd,line) ;
                        // on compte les espaces pour savoir combien on a de parametres
                        int nspace=0;
                        for(unsigned int j=0; j!=line.size(); ++j)
                            nspace+=( line.at(j)==' ');
                        int ntmp = nspace+1 ;
                        if (ntmp==1)
                        {
                            int* tmp = new int[ntmp] ;
                            // on met les entiers dans les cases
                            size_t spos ;
                            size_t isize ;
                            size_t sdeb = 0 ;
                            int ind = 0 ;
                            do {
                                spos = line.find(' ', sdeb) ;
                                isize = spos - sdeb ;
                                istringstream(line.substr(sdeb, isize)) >> tmp[ind] ;
                                ind ++;
                                sdeb = spos+1 ;
                            } while (spos!=string::npos ) ;
                            // on a toutes les composantes dans tmp
                            J[i]=tmp[0];
                        }
                    }
                }

                if( line.compare("$CData")==0 )
                {
                    fd>>data_size;
                    C_data.resize(data_size);
                    getline(fd,line);
                    for( i=0; i<data_size ; i++)
                    {
                        // on ne sais pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
                        getline(fd,line) ;
                        // on compte les espaces pour savoir combien on a de parametres
                        int nspace=0;
                        for(unsigned int j=0; j!=line.size(); ++j)
                            nspace+=( line.at(j)==' ');
                        int ntmp = nspace+1 ;
                        if (ntmp==6)
                        {
                            int* tmp = new int[ntmp] ;
                            // on met les entiers dans les cases
                            size_t spos ;
                            size_t isize ;
                            size_t sdeb = 0 ;
                            int ind = 0 ;
                            do {
                                spos = line.find(' ', sdeb) ;
                                isize = spos - sdeb ;
                                istringstream(line.substr(sdeb, isize)) >> tmp[ind] ;
                                ind ++;
                                sdeb = spos+1 ;
                            } while (spos!=string::npos ) ;
                            // on a toutes les composantes dans tmp
                            C_data.X[0][i]=C(tmp[0],tmp[3]);
                            C_data.X[1][i]=C(tmp[1],tmp[4]);
                            C_data.X[2][i]=C(tmp[2],tmp[5]);
                            //cout<<C_data.X[0][i]<<"  -   "<<C_data.X[1][i]<<"   -   "<<C_data.X[2][i]<<endl;
                        }
                    }
                }

            }
        }
    }
}

int MatSparseC3::EditProfileFromMesh(Maillage3D _maillage)
{
    size=_maillage.GetNodesSize();
    P.resize(size);

    map<int,set<int>*>::iterator map_it;
    set<int>::iterator set_it;
    int idx,i,j;
    int data_size =0;
    int data_cursor =0;
    int row_cursor=0;

    for(i=0;i<size;i++)
        P[i]=0;

    set<int>* neighbors_set;
    for (map_it =_maillage.GetNeighboursMap()->begin(); map_it !=_maillage.GetNeighboursMap()->end();++map_it)
    {
        idx=map_it->first;
        neighbors_set=map_it->second;
        if(neighbors_set!=0 && neighbors_set->size()>0)
        {
            for (set_it=neighbors_set->begin();set_it!=neighbors_set->end();++set_it)
                data_size++;
        }
    }

    for(j=0;j<3;j++)
    {
        C_data.X[j].resize(data_size);
    }
    J.resize(data_size);
    cout<<data_size<<" donnees non nulles dans le profil de la Matrice systeme "<<endl<<"Systeme a "<<100-((data_size*100.)/(size*size))<<"% sparse."<<endl;

    cout<<"Editing profile data of the Sparse Matrix to match the mesh configuration."<<endl;
    for (map_it =_maillage.GetNeighboursMap()->begin(); map_it !=_maillage.GetNeighboursMap()->end();++map_it)
    {
        idx=map_it->first;
        neighbors_set=map_it->second;
        if(neighbors_set!=0 && neighbors_set->size()>0)
        {
            for (set_it=neighbors_set->begin();set_it!=neighbors_set->end();++set_it)
            {
                J[data_cursor]=*set_it;
                for(i=row_cursor;i<size;i++)
                    P[i+1]++;
                data_cursor++;
            }
            row_cursor++;
        }
    }
    return 1;
}

int  MatSparseC3::IdxSearch(int _row,int _column)
{
    int j;
    //Search on the interval between the start of the current row and the next one
    // we use row-1 as starting index since the standard algebraic row numerotation starts at 1 while P[] is a "standard" array
    for(j=P[_row-1];j<min(long(P[_row]),J.size());j++)
    {
        if (J[j]==_column)
            return j;
    }
    return -1;
}

int MatSparseC3::MatDiagMul(MatDiag W)
{
    if(W.size()==size)
    {
        int i,j,k;
        for(k=0;k<3;k++)
        {
            for (i=0;i<size;i++)
            {
                for(j=P[i];j<P[i+1];j++)
                {

                    (C_data.X[k])[j]=(C_data.X[k])[j]*(W.X[j])[J[j]];
                }
            }
        }
    }
}

int MatSparseC3::DiagMatMul(MatDiag W)
{
    if(W.size()==size)
    {
        int i,k,j;
        for (j=0;j<3;j++)
        {
            for (i=0;i<size;i++)
            {
                for(k=P[i];k<P[i+1];k++)
                {

                    (C_data.X[j])[k]=(C_data.X[j])[k]*(W.X[j])[i];
                }
            }
        }
    }
}
void MatSparseC3::MatVectMul(VectorC3 & _X, VectorC3 &_o,int _comp)
{
    if(_X.size()==size)
    {
        int i,k,j;
        if(_comp ==-1 || _comp>2)
        {
            for (j=0;j<3;j++)
            {
                (_o.X[j]).resize(size);
                for (i=0;i<size;i++)
                {
                    (_o.X[j])[i]=(0,0);
                    for(k=P[i];k<min(long(P[i+1]),J.size());k++)
                    {
                        //Debug for insertions
                        //cout<<"Insert in ["<<i<<"] : "<<C((C_data.X[j])[k]*(_X.X[j])[J[k]])<< "at k = "<<k<<endl;
                        (_o.X[j])[i]+=C((C_data.X[j])[k]*(_X.X[j])[J[k]-1]);
                    }
                }
            }
        }
        else
        {
            j=_comp;
            (_o.X[j]).resize(size);
            for (i=0;i<size;i++)
            {
                (_o.X[j])[i]=(0,0);
                for(k=P[i];k<min(long(P[i+1]),J.size());k++)
                {
                    _o.X[j][i]+=C((C_data.X[j])[k]*(_X.X[j])[J[k]-1]);
                }
            }
        }
    }
}

