#include<complex>
#include<tuple>

#include "rhs.h"
#include "C3.h"


RHS::RHS(ifstream& fd, map<int, R3> *node_data)
{
    if (node_data!=0)
    {
        rhs = new vector<C3>;
        int size;
        vector<C>* rhsX= new vector<C>;
        vector<C>* rhsY= new vector<C>;
        vector<C>* rhsZ= new vector<C>;


        vector<double>* posX = new vector<double>;
        vector<double>* posY = new vector<double>;
        vector<double>* posZ = new vector<double>;
        string line;
        if(fd.is_open() ) {
            // on lit les lignes une par une et on avise en fonction des mots cle
            while ( fd.good() ) {
                getline (fd,line);
                // on regarde si on a un mot cle
                if( line[0] == '$' ) {

                    if( line.compare("$RHSExpression")==0 )
                    {
                        getline (fd,line);
                        if( line.compare("$XExpression")==0 )
                        {
                            getline(fd,line);
                            //save the data of the expression
                            expr[0]=line;
                        }
                        getline(fd,line);
                        if( line.compare("$YExpression")==0 )
                        {
                            getline(fd,line);
                            expr[1]=line;
                        }
                        getline(fd,line);
                        if( line.compare("$ZExpression")==0 )
                        {
                            getline(fd,line);
                            expr[2]=line;
                        }
                    }
                    else if( line.compare("$RHSExplicit")==0 )
                        //The node <-> Value is explicitely given
                    {
                        getline(fd,line);//Go to next line : the data corresponds to the size of the vector
                        fd>>size;
                        if (size==node_data->size()) //Input Data size corresponds to the Mesh Data size
                        {
                            while(line.size()!=0)
                            {
                                // on ne sait pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
                                getline(fd,line) ;
                                // on compte les espaces pour savoir combien on a de parametres
                                int nspace=0;
                                for(unsigned int j=0; j!=line.size(); ++j)
                                    nspace+=( line.at(j)==' ');
                                int ntmp = nspace+1 ;
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
                                // on a tous les entiers dans tmp
                                // Check if ntmp=7 (n°noeud + six composantes)
                                if(ntmp==7)
                                {
                                    //TODO
                                }
                                else
                                {
                                    cout<<"Discrepancy in data size. Abort"<<endl;
                                }
                            }
                        }
                        else
                        {
                            cout<<"Incompatible data size. Aborting."<<endl;
                            break;
                        }
                    }
                    else if( line.compare("$RHSImplicit")==0 )
                        //The given data only contains values, we'll implicitely assume the data is ordered according to the mesh's node ordering
                    {
                        while(line.size()!=0)
                        {
                            // on ne sait pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
                            getline(fd,line) ;
                            // on compte les espaces pour savoir combien on a de parametres
                            int nspace=0;
                            for(unsigned int j=0; j!=line.size(); ++j)
                                nspace+=( line.at(j)==' ');
                            int ntmp = nspace+1 ;
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
                            // on a tous les entiers dans tmp
                            // Check if ntmp=6 (six composantes)
                            if(ntmp==6)
                            {
                                //TODO
                            }
                            else
                            {
                                cout<<"Discrepancy in data size. Abort"<<endl;
                            }
                        }
                    }
                }
            }
        }
        delete rhsX;
        delete rhsY;
        delete rhsZ;
        delete posX;
        delete posY;
        delete posZ;
    }
}
