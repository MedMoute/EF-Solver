#include <QFileDialog>
#include <chrono>

#include "utils.h"
#include "map_mainwindow.h"
#include "ui_map_mainwindow.h"
#include "C3.h"

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
}

MAP_MainWindow::~MAP_MainWindow()
{
    delete ui;
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
        if (bc->GetSize()>0 || bc->GetExpr()[0].size()>0)
        {
            ui->radioButton_homogeneous->setChecked(false);
        }
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
        ui->radioButton_homogeneous->setChecked(true);
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

    problem = new FEProblem(maillage,rhs,bc);

    ui->pushButton_solve->setEnabled(true);
    ui->checkBox_dir->setEnabled(true);
    ui->checkBox_it->setEnabled(true);
    ui->comboBox_solver->setEnabled(true);
}

void MAP_MainWindow::on_pushButton_solve_released()
{
    util::print_separator();
    cout<<"Lancement de la routine de resolution du probleme."<<endl;
    solver = new Solver(problem,GC,0.0000001);
    if (solver->GetRes()==1)
        ui->pushButton_save_output->setEnabled(true);
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

    if (solver->GetRes()==1)
        ui->pushButton_save_output->setEnabled(true);
}

void MAP_MainWindow::on_pushButton_save_output_released()
{
//Open dialog in order to save the Output vector

}
