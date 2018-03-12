#ifndef MAP_MAINWINDOW_H
#define MAP_MAINWINDOW_H

#include <QWidget>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "maillage3d.h"
#include "rhs.h"
#include "bc.h"
#include "feproblem.h"
#include "solver.h"

namespace Ui {
class MAP_MainWindow;
}

class MAP_MainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit MAP_MainWindow(QWidget *parent = 0);

    void AssembleAll();
    ~MAP_MainWindow();

public slots:
    void on_pushButton_load_3D_released();
    void on_pushButton_load_RHS_released();
    void on_pushButton_mat_assemb_released();
    void on_pushButton_load_BC_released();
    void on_pushButton_solve_released();
    void on_pushButton_load_Mat_released();

private:
    Ui::MAP_MainWindow *ui;
    Maillage3D *maillage;
    RHS *rhs;
    BC *bc;
    FEProblem *problem;
    Solver *solver;
};

#endif // MAP_MAINWINDOW_H
