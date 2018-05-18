#include "parametersdialog.h"
#include "ui_parametersdialog.h"
#include <algorithm>
ParametersDialog::ParametersDialog(std::vector<int>*partitions, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ParametersDialog)
{
    ui->setupUi(this);
    for (unsigned int i=0;i<partitions->size();i++)
    {
        ui->comboBox->addItem(QString::fromStdString(std::to_string(partitions->at(i))));
    }
    ui->comboBox->setEnabled(false);
}

void ParametersDialog::on_pushButton_validate_released()
{
    if(!ui->radioButton_global->isChecked())//Partitionned parametrisation
    {
        part_data.insert(std::pair<int,double>(ui->comboBox->currentText().toInt(),ui->doubleSpinBox->value()));
    }
}

void ParametersDialog::on_radioButton_global_toggled()
{
    if(ui->radioButton_global->isChecked())
        ui->comboBox->setEnabled(false);
    else
        ui->comboBox->setEnabled(true);
}

ParametersDialog::~ParametersDialog()
{
    delete ui;
}
