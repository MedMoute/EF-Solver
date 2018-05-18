#include <chrono>

#include <QFileDialog>
#include <QDir>
#include <QtWidgets/QMessageBox>
#include <QtGui/QScreen>
#include <QHBoxLayout>
#include <QVBoxLayout>

#include "map_mainwindow.h"
#include "ui_map_mainwindow.h"
#include "C3.h"

#include "viewer3d.h"
#include "scatterdatamodifier.h"

MAP_MainWindow::MAP_MainWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MAP_MainWindow)
{
    ui->setupUi(this);
    //INIT member pointers to NULL
    rhs=0;
    bc=0;
    maillage=0;
    problem=0;
    solver=0;
    m_meshVector=new vector<Triangle>;
}

MAP_MainWindow::~MAP_MainWindow()
{
    delete ui;
    delete m_meshVector;
}

void MAP_MainWindow::on_pushButton_load_BC_released()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "",
                                                    tr("Array Data File (*.arr);;Any files(*)"));
    std::ifstream FILE;
    FILE.open(fileName.toStdString().c_str(), std::ios::in);

    if (FILE.fail())
    {
        std::cout<<"Erreur lors de l'ouverture de FILE"<<std::endl;
    }

    else{
        /* Creation des donnees de conditions au bord a partir du fichier lu */
        util::print_separator();

        auto start_time = std::chrono::system_clock::now();
        bc =new BC(FILE,maillage->GetBoundaryNodesMap());
        auto end_time = std::chrono::system_clock::now();
        auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
        std::cout<<"Chargement du fichier des conditions aux bord : "<<fileName.toStdString()<<" effectue en "<<elapsed_seconds<<"s."<<std::endl;
    }
}

void MAP_MainWindow::on_pushButton_load_3D_released()
{
    //reset RHS and BC data
    if (rhs!=0)
    {
        delete rhs;
    }
    if (bc!=0)
    {
        delete bc;
    }
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "",
                                                    tr("GMSH Mesh File (*.msh);;Any files(*)"));
    std::ifstream FILE;
    FILE.open(fileName.toStdString().c_str(), std::ios::in);

    if (FILE.fail())
    {
        std::cout<<"Erreur lors de l'ouverture de FILE"<<std::endl;
    }

    else{
        /* Creation des donnees de maillage a partir du fichier lu */

        auto start_time = std::chrono::system_clock::now();
        maillage =new Maillage3D(FILE);
        auto end_time = std::chrono::system_clock::now();
        auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
        std::cout<<"Chargement du fichier "<<fileName.toStdString()<<" effectue en "<<elapsed_seconds<<"s."<<std::endl;

        maillage->DisplayBasicData();
        maillage->DisplayBoundaryElements();

        ui->comboBox_problem->setEnabled(true);
        ui->pushButton_mat_assemb->setEnabled(true);
        ui->pushButton_load_BC->setEnabled(true);
        ui->pushButton_load_RHS->setEnabled(true);
        ui->pushButton_vis_3D->setEnabled(true);
    }
}

void MAP_MainWindow::  on_pushButton_load_RHS_released()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "",
                                                    tr("Array Data File (*.arr);;Any files(*)"));
    std::ifstream FILE;
    FILE.open(fileName.toStdString().c_str(), std::ios::in);

    if (FILE.fail())
    {
        std::cout<<"Erreur lors de l'ouverture de FILE"<<std::endl;
    }

    else{
        /* Creation des donnees du terme de second membre a partir du fichier lu */
        util::print_separator();

        auto start_time = std::chrono::system_clock::now();
        rhs =new RHS(FILE,maillage->GetNodesMap());
        auto end_time = std::chrono::system_clock::now();
        auto elapsed_seconds = std::chrono::duration<float>(end_time-start_time).count();
        std::cout<<"Chargement du fichier du second terme : "<<fileName.toStdString()<<" effectue en "<<elapsed_seconds<<"s."<<std::endl;

    }
}

void MAP_MainWindow::on_pushButton_mat_assemb_released()
{
    if (problem)
        delete problem;

    util::print_separator();
    cout<<"Lancement de la routine d'assemblage des Matrices"<<endl;
    //map<int,tuple<int,int,int,int>>* m_ptr = maillage->GetElementsMap();
    //DisplayElementaryMassMatrix(m_ptr->begin()->first);
    //DisplayElementaryStiffnessMatrix(m_ptr->begin()->first);
    int op=ui->comboBox_problem->currentIndex();

    problem = new FEProblem(maillage,rhs,bc,op,ui->checkBox_unitary->isChecked());

    ui->pushButton_solve->setEnabled(true);
    ui->checkBox_dir->setEnabled(true);
    ui->checkBox_it->setEnabled(true);
    ui->comboBox_solver->setEnabled(true);
}

void MAP_MainWindow::on_pushButton_solve_released()
{

    util::print_separator();
    cout<<"Lancement de la routine de resolution du probleme."<<endl;
    solver = new Solver(problem,SolverMethod(ui->comboBox_solver->currentIndex()+1),0.0000001);

    ui->pushButton_save_output->setEnabled(true);

cout<<"Erreur finale : "<<solver->computeError()<<endl;
}

void MAP_MainWindow::on_pushButton_load_Mat_released()
{
    util::print_separator();
    cout<<"Chargement direct d'une matrice Sparse"<<endl;
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Matrix File"),
                                                    "",
                                                    tr("Sparse Matrix Data File (*.spa);;Any files(*)"));
    std::ifstream FILE;
    FILE.open(fileName.toStdString().c_str(), std::ios::in);

    if (FILE.fail())
    {
        std::cout<<"Erreur lors de l'ouverture de FILE"<<std::endl;
    }
    MatSparseC3 mat(FILE);
    FILE.close();

    cout<<"Chargement direct du second membre du systeme lineaire"<<endl;
    fileName = QFileDialog::getOpenFileName(this, tr("Open Vector File"),
                                            "",
                                            tr("Vector Data File (*.vec);;Any files(*)"));
    FILE.open(fileName.toStdString().c_str(), std::ios::in);

    if (FILE.fail())
    {
        std::cout<<"Erreur lors de l'ouverture de FILE"<<std::endl;
    }

    VectorC3 vec(FILE);
    FILE.close();
    ui->pushButton_solve->setEnabled(true);
    ui->checkBox_dir->setEnabled(true);
    ui->checkBox_it->setEnabled(true);
    ui->comboBox_solver->setEnabled(true);

    solver = new Solver(mat,vec,SolverMethod(ui->comboBox_solver_2->currentIndex()+1),0.00001);
    cout<<"Erreur finale : "<<solver->computeError(mat,vec)<<endl;
    solver->displaySolution();
}

void MAP_MainWindow::on_pushButton_save_output_released()
{
    //Open dialog in order to save the Output vector

}

void MAP_MainWindow::on_pushButton_vis_3D_released()
{
    if (maillage)
    {
        string path = (QDir::currentPath()+"tmp.obj").toStdString();
        int res = maillage->ExportMeshAsOBJFile(path);
        if (res==0)
        {
            Viewer3D* m_viewer3D = new Viewer3D(this);

            CreateMeshBuffer();

            m_viewer3D->sceneModifier()->addTriangleMeshCustomMaterial(QString::fromStdString(path),m_meshVector);
            m_viewer3D->show();
        }
    }
}

void MAP_MainWindow::CreateMeshBuffer()
{
    map<int,std::tuple<int,int,int>>* face_m = maillage->GetFacesMap();
    map<int,std::tuple<int,int,int>>::iterator face_it;
    map<int,R3>* node_m = maillage->GetNodesMap();
    map<int,R3>::iterator node_it;

    for(face_it=face_m->begin();face_it!=face_m->end();++face_it)
    {
        Triangle tmp_tri;
        int N1=std::get<0>((*face_it).second);
        node_it =node_m->find(N1);
        tmp_tri.vertices[0].p.setX((*node_it).second.X_());
        tmp_tri.vertices[0].p.setY((*node_it).second.Y_());
        tmp_tri.vertices[0].p.setZ((*node_it).second.Z_());
        int N2=std::get<1>((*face_it).second);
        node_it =node_m->find(N2);
        tmp_tri.vertices[1].p.setX((*node_it).second.X_());
        tmp_tri.vertices[1].p.setY((*node_it).second.Y_());
        tmp_tri.vertices[1].p.setZ((*node_it).second.Z_());
        int N3=std::get<2>((*face_it).second);
        node_it =node_m->find(N3);
        tmp_tri.vertices[2].p.setX((*node_it).second.X_());
        tmp_tri.vertices[2].p.setY((*node_it).second.Y_());
        tmp_tri.vertices[2].p.setZ((*node_it).second.Z_());
        m_meshVector->push_back(tmp_tri);
    }
}

