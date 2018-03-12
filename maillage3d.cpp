#include "maillage3D.h"
#include "utils.h"

#include <chrono>



Maillage3D::Maillage3D(ifstream& fd)
{
    nodes = new map<int,R3> ;
    boundary_nodes = new map<int,R3>;
    partitions = new vector<int>;
    edges = new map<int,pair<int,int>> ;
    faces = new map<int,tuple<int,int,int>> ;
    faces_bord = new map<int,tuple<int,int,int>>;
    elements = new map<int,tuple<int,int,int,int>>;
    elements_partition = new map<int,int> ;

    neighbors = new map<int,set<int>*>;
    neighbors_faces = new map<int,set<int>*>;
    neighbors_elements = new map<int,set<int>*>;

    boundary_neighbors = new map<int,set<int>*>;
    neighbors_interface = new map<int,pair<int,int>>;
    faces_bord_partitions = new map<int,int>;
    triangles_interfaces =new map<int,pair<int,int>> ;
    triangles_interfaces_nodes = new map<int,tuple<int,int,int>> ;

    // POUR LES FICHIERS VERSION 2.2+ SEULEMENT
    string line; //dump is a string to empty the line data after we've extracted the data from line.
    // on lit dans le fichier fd
    if(fd.is_open() ) {
        // on lit les lignes une par une et on avise en fonction des mots cle
        while ( fd.good() ) {
            getline (fd,line);
            // on regarde si on a un mot cle
            if( line[0] == '$' ) {

                if( line.compare("$MeshFormat")==0 ) {
                    // je lis le format mais je n'en fais rien
                    getline(fd,line) ;
                    // je lis le mot cle de fin
                    getline(fd,line) ;
                    if( ! line.compare("$EndMeshFormat")==0 ) {
                        cout << "Probleme de lecture du fichier maillage" << endl ;
                        cout << "Mot cle MeshFormat non suivi de EndMeshFormat" << endl ;
                        abort() ;
                    }
                }
                else if( line.compare("$Nodes")==0 ) {
                    // lecture du nombre de noeuds
                    fd >> n_nodes;

                    // on lit les coords
                    for (int i = 0; i< n_nodes ; i++)
                    {
                        // variable dummy pour le numero de noeud
                        int idx;
                        float x,y,z;

                        fd >> idx >> x >> y>>z ;
                        fd.ignore(numeric_limits<streamsize>::max(), '\n');//ignore data until next line break
                        nodes->insert(vertexR3(idx,R3(x,y,z)));
                        if(i%10000==0)
                        {
                            cout<<"["<<int(i)*100/n_nodes<<"%]"<<"Analysing node #"<<i<<"/"<<n_nodes<<endl;
                        }
                    }
                    cout<<"Node Analysis finished"<<endl;
                    util::print_separator();

                    getline(fd,line) ;
                    // on lit le mot cle de fin

                    if( ! line.compare("$EndNodes")==0 ) {
                        cout << "Probleme de lecture du fichier maillage" << endl ;
                        cout << "Mot cle Nodes non suivi de EndNodes" << endl ;
                        abort() ;
                    }

                }
                else if( line.compare("$Elements")==0 ) {
                    // on lit le nb d'elements
                    fd >> n_elems;
                    // on passe a la ligne pour la suite
                    getline(fd,line) ;

                    for( int i=0; i<n_elems ; i++)
                    {


                        // on ne sais pas combien il y a de parametres sur chaque ligne : on lit ligne par ligne
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
                        // on remplit les membres de maillage
                        int nb_sommet(util::type2nnodes(tmp[1]));  // nb de sommets pour cet element

                        if (nb_sommet==2)//l'element est un segment
                        {
                            edges->insert(edgeR3(tmp[0],pair<int,int>(tmp[ntmp-2],tmp[ntmp-1])));
                        }
                        else if(nb_sommet==3)//l'element est un triangle
                        {
                            faces->insert(faceR3(tmp[0],tuple<int,int,int>(tmp[ntmp-3],tmp[ntmp-2],tmp[ntmp-1])));

                            //Remplissage de faces_bord_partitions
                            int j;
                            int n_interf_proches=0;
                            for (j=0;j<ntmp;j++)
                            {
                                if (tmp[j]<0)
                                {
                                    n_interf_proches++;
                                }
                            }
                            if (n_interf_proches ==0)//L'element traité n'est pas sur une interface, il est donc sur le bord du domaine
                            {
                                faces_bord_partitions->insert(pair<int,int>(tmp[0],tmp[ntmp-4]));
                                faces_bord->insert(faceR3(tmp[0],tuple<int,int,int>(tmp[ntmp-3],tmp[ntmp-2],tmp[ntmp-1])));
                                //On profite des information actuellement disponibles pour remplir la map des boundary_nodes.
                                int k;
                                for (k=1;k<=3;k++)
                                {
                                    map<int,R3>::iterator nodes_it=nodes->find(tmp[ntmp-k]);
                                    if(nodes_it!=nodes->end())
                                    {
                                        boundary_nodes->insert(pair<int,R3>(nodes_it->first,nodes_it->second));
                                    }
                                }
                            }
                            //remplissage neighbor elements
                            for(j=1;j<=3;j++)
                            {
                                int idx = tmp[ntmp-j];
                                map<int,set<int>*>::iterator it_sets;
                                it_sets =neighbors_faces->find(idx);
                                set<int>* idx_neigh;
                                if (it_sets!=neighbors_faces->end())//check if the neighor set already exists

                                {//If yes fetch it
                                    idx_neigh = it_sets->second;
                                }
                                else
                                {//Else add the previously created set to the map of neigbours
                                    idx_neigh=new set<int>;
                                    neighbors_faces->insert(pair<int,set<int>*>(idx,idx_neigh));
                                }
                                //Then add the element to the set corresponding to the node
                                idx_neigh->insert(tmp[0]);
                            }

                        }
                        else if(nb_sommet==4)//l'element est un tetraedre
                        {
                            elements->insert(elementR3(tmp[0],tuple<int,int,int,int>(tmp[ntmp-4],tmp[ntmp-3],tmp[ntmp-2],tmp[ntmp-1])));

                            //remplissage neighbor elements
                            int j;
                            for(j=1;j<=4;j++)
                            {
                                int idx = tmp[ntmp-j];
                                map<int,set<int>*>::iterator it_sets;
                                it_sets =neighbors_elements->find(idx);
                                set<int>* idx_neigh;
                                if (it_sets!=neighbors_elements->end())//check if the neighor set already exists

                                {//If yes fetch it
                                    idx_neigh = neighbors_elements->find(idx)->second;
                                }
                                else
                                {//Else add the previously created set to the map of neigbours
                                    idx_neigh=new set<int>;
                                    neighbors_elements->insert(pair<int,set<int>*>(idx,idx_neigh));
                                }
                                //Then add the element to the set corresponding to the node
                                idx_neigh->insert(tmp[0]);
                            }

                            //Remplissage de voisins interface
                            // Supposition effectuée:
                            // Les données négatives dans tmp correspondent aux partitions
                            // dont l'element est à proximité, qui ne SONT PAS la partition de l'element.
                            int n_part_proches=0;
                            for (j=0;j<ntmp;j++)
                            {
                                if (tmp[j]<0)
                                {
                                    n_part_proches++;
                                }
                            }

                            //Ajout de l'information de la partition de l'element
                            elements_partition->insert(pair<int,int>(tmp[0],tmp[ntmp-nb_sommet-n_part_proches-1]));

                            //Insertions dans neigbors_interface
                            for (j=0;j<ntmp;j++)
                            {
                                if (tmp[j]<0)
                                {
                                    neighbors_interface->insert(pair<int,pair<int,int>>(tmp[0],pair<int,int>(tmp[ntmp-nb_sommet-n_part_proches-1],-1*tmp[j])));
                                }
                            }
                            //on ajoute la partition a la liste
                            partitions->push_back(tmp[ntmp-nb_sommet-n_part_proches-1]);

                        }

                        delete [] tmp ;
                        if(i%10000==0)
                        {
                            cout<<"["<<int(i)*100/n_elems<<"%]"<<"Analysing element #"<<i<<"/"<<n_elems<<endl;
                        }
                    }
                    cout<<"Element Analysis finished"<<endl;
                    util::print_separator();

                    getline(fd,line) ;
                    if( ! line.compare("$EndElements")==0 ) {
                        cout << "Probleme de lecture du fichier maillage" << endl ;
                        cout << "Mot cle Nodes non suivi de EndElements" << endl ;
                    }

                }
            }
        }
    }
    //Traitement des voisinages de noeuds
    map<int,tuple<int,int,int,int>>::iterator it;
    for (it=elements->begin(); it!=elements->end(); ++it)
    {
        int idx1 = get<0>(it->second);
        int idx2 = get<1>(it->second);
        int idx3 = get<2>(it->second);
        int idx4 = get<3>(it->second);

        map<int,set<int>*>::iterator it_sets;
        it_sets =neighbors->find(idx1);
        set<int>* idx1_neigh;
        if (it_sets!=neighbors->end())//check if the neighor set already exists

        {//If yes fetch it
            idx1_neigh = neighbors->find(idx1)->second;
        }
        else
        {//Else add the previously created set to the map of neigbours
            idx1_neigh=new set<int>;
            neighbors->insert(pair<int,set<int>*>(idx1,idx1_neigh));
        }
        //Then add the neighbors to the set
        idx1_neigh->insert(idx2);
        idx1_neigh->insert(idx3);
        idx1_neigh->insert(idx4);
        //Adding the node itself for it to be considered in the Matrix Construction
        idx1_neigh->insert(idx1);



        //Do the same for the other nodes of the element

        it_sets =neighbors->find(idx2);
        set<int>* idx2_neigh;
        if (it_sets!=neighbors->end())
        {idx2_neigh = neighbors->find(idx2)->second;}
        else
        {idx2_neigh= new set<int>;neighbors->insert(pair<int,set<int>*>(idx2,idx2_neigh));}
        idx2_neigh->insert(idx3);
        idx2_neigh->insert(idx1);
        idx2_neigh->insert(idx4);
        idx2_neigh->insert(idx2);

        it_sets =neighbors->find(idx3);
        set<int>* idx3_neigh;
        if (it_sets!=neighbors->end())
        {idx3_neigh = neighbors->find(idx3)->second;}
        else
        {idx3_neigh=new set<int>;neighbors->insert(pair<int,set<int>*>(idx3,idx3_neigh));}
        idx3_neigh->insert(idx2);
        idx3_neigh->insert(idx1);
        idx3_neigh->insert(idx4);
        idx3_neigh->insert(idx3);


        it_sets =neighbors->find(idx4);
        set<int>* idx4_neigh;
        if (it_sets!=neighbors->end())
        {idx4_neigh = neighbors->find(idx4)->second;}
        else
        {idx4_neigh=new set<int>;neighbors->insert(pair<int,set<int>*>(idx4,idx4_neigh));}
        idx4_neigh->insert(idx1);
        idx4_neigh->insert(idx2);
        idx4_neigh->insert(idx3);
        idx4_neigh->insert(idx4);

    }

    //Traitement des voisinages le long des triangles du bord du domaine
    map<int,tuple<int,int,int>>::iterator it_b;
    for (it_b=faces_bord->begin(); it_b!=faces_bord->end(); ++it_b)
    {
        int idx1 = get<0>(it_b->second);
        int idx2 = get<1>(it_b->second);
        int idx3 = get<2>(it_b->second);

        map<int,set<int>*>::iterator it_sets;
        it_sets =boundary_neighbors->find(idx1);
        set<int>* idx1_neigh;
        if (it_sets!=boundary_neighbors->end())//check if the neighor set already exists

        {//If yes fetch it
            idx1_neigh = boundary_neighbors->find(idx1)->second;
        }
        else
        {//Else add the previously created set to the map of neigbours
            idx1_neigh=new set<int>;
            boundary_neighbors->insert(pair<int,set<int>*>(idx1,idx1_neigh));
        }
        //Then add the neighbors to the set
        idx1_neigh->insert(idx2);
        idx1_neigh->insert(idx3);

        //Do the same for the other nodes of the element

        it_sets =boundary_neighbors->find(idx2);
        set<int>* idx2_neigh;
        if (it_sets!=boundary_neighbors->end())
        {idx2_neigh = boundary_neighbors->find(idx2)->second;}
        else
        {idx2_neigh = new set<int>;boundary_neighbors->insert(pair<int,set<int>*>(idx2,idx2_neigh));}
        idx2_neigh->insert(idx3);
        idx2_neigh->insert(idx1);

        it_sets =boundary_neighbors->find(idx3);
        set<int>* idx3_neigh=new set<int>;
        if (it_sets!=boundary_neighbors->end())
        {idx3_neigh = boundary_neighbors->find(idx3)->second;}
        else
        {idx3_neigh=new set<int>;boundary_neighbors->insert(pair<int,set<int>*>(idx3,idx3_neigh));}
        idx3_neigh->insert(idx2);
        idx3_neigh->insert(idx1);
    }

    //Traitement du vecteur partitions (suppresion des doublons)

    vector< int >::iterator r , w ;
    set< int > tmpset ;

    for( r = partitions->begin() , w = partitions->begin() ; r != partitions->end() ; ++r )
    {
        if( tmpset.insert( *r ).second )
        {
            *w++ = *r ;
        }
    }
    partitions->erase( w , partitions->end() );
    //On ajoute le traitement des interfaces et la création des structures de triangles posés sur les interfaces
    if (partitions->size()>1)
        CreateInterfaceTriangles();
}

int Maillage3D::GetPartitionSize( int _part)
{
    //TODO
    int res=0;
    return res;
}


void Maillage3D::CreateInterfaceTriangles()
{
    cout<<"Creating Interfaces between partitions. This step can take a while..."<<endl;

    int i,j;
    int part_i,part_j;
    int elem_i,elem_j;
    tuple<int,int,int,int> nodes_i,nodes_j;
    map<int,tuple<int,int,int,int>>::iterator it_nodes_i,it_nodes_j;
    //Debug for CPU time
    auto start_time = chrono::system_clock::now();
    auto elapsed_seconds_common = chrono::duration<float>(start_time-start_time).count();
    int commonTriCount=0;
    int neigh_i_size=0;
    int neigh_j_size=0;

    for (i=0;i<partitions->size();i++)
    {
        part_i=partitions->at(i);
        for (j=0;j<i;j++)
        {
            part_j=partitions->at(j);
            if (part_i!=part_j)
            {
                map<int,pair<int,int>>::iterator it_i;
                map<int,pair<int,int>>::iterator it_j;

                //Parcours des éléments à proximité d'une interface et recherche de l'interface en cours
                for (it_i = neighbors_interface->begin(); it_i != neighbors_interface->end(); ++it_i )
                    if (it_i->second == pair<int,int>(part_i,part_j))
                    {//Parcours des elements de la partition i à proximité de j
                        elem_i=it_i->first;
                        //recuperation des noeuds de l'elem de i en cours de traitement
                        it_nodes_i=elements->find(elem_i);
                        neigh_i_size++;
                        if (it_nodes_i!=elements->end())
                        {
                            nodes_i=it_nodes_i->second;
                        }
                        for (it_j = neighbors_interface->begin(); it_j != neighbors_interface->end(); ++it_j )
                            if (it_j->second == pair<int,int>(part_j,part_i))
                            {//Parcours des elements de la partition j à proximité de i
                                elem_j=it_j->first;
                                //recuperation des noeuds de l'elem de j en cours de traitemetn
                                it_nodes_j=elements->find(elem_j);
                                neigh_j_size++;
                                if (it_nodes_j!=elements->end())
                                {
                                    nodes_j=it_nodes_j->second;
                                    pair<bool,tuple<int,int,int>> Tri;
                                    auto start_time_common = chrono::system_clock::now();
                                    Tri = CommonTriangle(nodes_i,nodes_j);
                                    auto end_time_common = chrono::system_clock::now();
                                    elapsed_seconds_common += chrono::duration<float>(end_time_common-start_time_common).count();
                                    commonTriCount++;
                                    if(Tri.first==1)
                                    {
                                        //cout<<get<0>(Tri.second)<<"-"<<get<1>(Tri.second)<<"-"<<get<2>(Tri.second)<<endl;
                                        triangles_interfaces->insert(pair<int,pair<int,int>>(triangles_interfaces->size()+1,pair<int,int>(part_i,part_j)));
                                        triangles_interfaces_nodes->insert(pair<int,tuple<int,int,int>>(triangles_interfaces_nodes->size()+1,tuple<int,int,int>(Tri.second)));
                                    }
                                }
                            }

                    }
            }
        }
        cout<<endl;
    }
    auto end_time = chrono::system_clock::now();
    auto elapsed_seconds = chrono::duration<float>(end_time-start_time).count();
    cout<<"Traitement et extractions des interfaces effectue en "<<elapsed_seconds<<"s."<<endl;
    cout<<commonTriCount<<" passes de commonTriangle effectuees en "<<elapsed_seconds_common<<"s."<<endl;
    cout<< "Duree moy. de commonTriangle: "<<elapsed_seconds_common/commonTriCount<<"s."<<endl;
    cout<<"Loop 1 size : "<< neigh_i_size<<" --- Loop 2 size : "<<neigh_j_size/neigh_i_size<<endl;
}

pair<bool,tuple<int,int,int>> Maillage3D::CommonTriangle(tuple<int,int,int,int> elem_i,tuple<int,int,int,int> elem_j)
{
                              tuple<int,int,int> a{-1,-1,-1};
int i_1=get<0>(elem_i);
int i_2=get<1>(elem_i);
int i_3=get<2>(elem_i);
int i_4=get<3>(elem_i);
int j_1=get<0>(elem_j);
int j_2=get<1>(elem_j);
int j_3=get<2>(elem_j);
int j_4=get<3>(elem_j);

vector<int> elem_i_s{i_1,i_2,i_3,i_4};
sort(elem_i_s.begin(),elem_i_s.end());
vector<int> elem_j_s{j_1,j_2,j_3,j_4};
sort(elem_j_s.begin(),elem_j_s.end());

int i_s_1=elem_i_s.at(0);
int i_s_2=elem_i_s.at(1);
int i_s_3=elem_i_s.at(2);
int i_s_4=elem_i_s.at(3);
int j_s_1=elem_j_s.at(0);
int j_s_2=elem_j_s.at(1);
int j_s_3=elem_j_s.at(2);
int j_s_4=elem_j_s.at(3);
//Debug output of commonTriangle
//cout<<"TEST "<<i_s_1<<"-"<<i_s_2<<"-"<<i_s_3<<"-"<<i_s_4<<" // " <<j_s_1<<"-"<<j_s_2<<"-"<<j_s_3<<"-"<<j_s_4<<endl;

//There are 14 permutation cases possible for 3 in 4 for sorted sets

if (i_s_1==j_s_1 && i_s_2==j_s_2 && i_s_3==j_s_3)// Cas |||o
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_2,i_s_3));
}
else if (i_s_1==j_s_1 && i_s_2==j_s_2 && i_s_4==j_s_4)// Cas ||o|
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_2,i_s_4));
}
else if (i_s_1==j_s_1 && i_s_3==j_s_3 && i_s_4==j_s_4)//cas |o||
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_3,i_s_4));
}
else if (i_s_2==j_s_2 && i_s_3==j_s_3 && i_s_4==j_s_4)// cas o|||
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_2,i_s_2,i_s_4));
}
else if (i_s_1==j_s_1 && i_s_2==j_s_2 && i_s_4==j_s_3)//cas  ||/
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_2,i_s_4));
}
else if (i_s_1==j_s_1 && i_s_2==j_s_2 && i_s_3==j_s_4)//cas  ||\
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_2,i_s_3));
}
else if (i_s_2==j_s_1 && i_s_3==j_s_3 && i_s_4==j_s_4)//cas  /||
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_2,i_s_3,i_s_4));
}
else if (i_s_1==j_s_2 && i_s_3==j_s_3 && i_s_4==j_s_4)//cas  \||
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_3,i_s_4));
}
else if (i_s_1==j_s_1 && i_s_2==j_s_3 && i_s_4==j_s_4)//cas |\|
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_2,i_s_4));
}
else if (i_s_1==j_s_1 && i_s_3==j_s_2 && i_s_4==j_s_4)// cas |/|
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_3,i_s_4));
}
else if (i_s_1==j_s_1 && i_s_3==j_s_2 && i_s_4==j_s_3)// cas |//
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_3,i_s_4));
}
else if (i_s_1==j_s_1 && i_s_2==j_s_3 && i_s_3==j_s_4)// cas | \ \
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_2,i_s_3));
}
else if (i_s_1==j_s_2 && i_s_2==j_s_3 && i_s_4==j_s_4)//cas \\|
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_1,i_s_2,i_s_4));
}
else if (i_s_2==j_s_1 && i_s_3==j_s_2 && i_s_4==j_s_4)//cas //|
//placeholder
{
    return pair<bool,tuple<int,int,int>>(true,tuple<int,int,int>(i_s_2,i_s_3,i_s_4));
}

else
{
return pair<bool,tuple<int,int,int>>(false,a);
}
}

void Maillage3D::DisplayBoundaryElements()
{
    util::print_separator();
    cout<<"Elements situes sur le bord du domaine : "<<endl;
    cout<<faces_bord_partitions->size()<<" elements sur le bord du maillage"<<endl;

}

void Maillage3D::DisplayInterfacesElements()
{
    util::print_separator();
    cout<<"Elements a proximite des interfaces : "<<endl;
    cout<<neighbors_interface->size()<<" elements proches d'une interface dans le maillage"<<endl;
    int i,j;
    int part_i,part_j;
    for (i=0;i<partitions->size();i++)
    {
        part_i=partitions->at(i);
        for (j=0;j<partitions->size();j++)
        {
            part_j=partitions->at(j);
            if (part_i!=part_j)
            {
                map<int,pair<int,int>>::iterator it;
                cout<<"Interface "<<part_i<<"/"<<part_j<<":"<<endl;
                //Parcours des éléments à proximité d'une interface et recherche de l'interface en cours
                for (it = neighbors_interface->begin(); it != neighbors_interface->end(); ++it )
                    if (it->second == pair<int,int>(part_i,part_j))
                    {
                        cout<<it->first<<" -";
                    }
            }
        }
        cout<<endl;
    }
    util::print_separator();
}

void Maillage3D::DisplayBasicData()
{

    util::print_separator();
    cout<<"Le maillage contient :"<<endl<<"* "<<this->GetNodesSize()<<"nodes"<<endl;
    cout<<"* "<<this->GetEdgesSize() <<" edges"<<endl;
    cout<<"* "<<this->GetFacesSize() <<" faces"<<endl;
    cout<<"* "<<this->GetElementsSize() <<" elements"<<endl;
    cout<<"On a trouve "<<triangles_interfaces->size()<<" triangles poses sur des interfaces"<<endl;

}

Maillage3D::~Maillage3D()
{}


