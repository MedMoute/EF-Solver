#ifndef PARAMETERSDIALOG_H
#define PARAMETERSDIALOG_H

#include <QDialog>
#include <vector>

namespace Ui {
class ParametersDialog;
}

class ParametersDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ParametersDialog(std::vector<int>* partitions,int _op,QWidget *parent = 0);
    ~ParametersDialog();
    Ui::ParametersDialog* GetUI(){return ui;}
    std::map<int,double> GetData(){return part_data;}
    int GetOP(){return op;}

private:
    Ui::ParametersDialog *ui;
    std::map<int,double> part_data;
    int op;

public slots:
    void on_pushButton_validate_released();
    void on_radioButton_global_toggled();

};

#endif // PARAMETERSDIALOG_H
