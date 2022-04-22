#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
#include <QLocale>
#include <locale.h>
#include <math.h>       /* round, floor, ceil, trunc */
MainWindow::MainWindow(QWidget *parent) :

    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    setlocale(LC_ALL, "C");
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_2_clicked()
{
    bool cali_status = ui->calibration_check->isChecked();
    int answer;
    if(cali_status){
        answer = QMessageBox::question(this,"","Finish this run?",QMessageBox::Yes,QMessageBox::No);
    }
    else{
        answer = QMessageBox::question(this,"","<b>Warning</b>: calibration might not exist.\nFinish run anyway?",QMessageBox::Yes,QMessageBox::No);
    }

    if(answer==QMessageBox::Yes){

        QMessageBox::about(this,"","New Run!");
        int runNo = std::stoi(ui->run->text().toStdString());
        runNo++;
        std::string newRun = std::to_string(runNo);

        ui->run->setText(QString::fromStdString(newRun));
        ui->subrun->setText("0");

        ui->calibration_check->setChecked(false);
    }

}

void MainWindow::on_button_movefile_clicked()
{
    QMessageBox::about(this,"","File was moved, new subrun");

    // This take the necessary info to create the folder and files
    int runNo = std::stoi(ui->run->text().toStdString());
    int subRunNo = std::stoi(ui->subrun->text().toStdString());
    double voltage = std::stod(ui->voltage->text().toStdString());
    double threshold = std::stod(ui->threshold->text().toStdString());
    std::string triggerCh = ui->trigger_channel->text().toStdString();
    std::string extra = ui->extra->text().toStdString();
    std::string primary = ui->primary_name->text().toStdString();
    //

    // here we create the folder and the move the files
    move_data_file(runNo,subRunNo,voltage,threshold,triggerCh,extra,primary);

    // update subrun number
    subRunNo++;
    std::string newSubRun = std::to_string(subRunNo);

    ui->subrun->setText(QString::fromStdString(newSubRun));
}

std::string MainWindow::folder_name(int run, int subrun, double voltage, double threshold, std::string triggerCh){
    std::string folder = "";
    std::string voltageS = changeVoltage(voltage);
    folder = "run" + std::to_string(run) + "_" + voltageS + "_" + std::to_string(static_cast<int>(threshold)) + "ADC_" + triggerCh + "/";
    return folder;
}

bool MainWindow::move_data_file(int run, int subrun, double voltage, double threshold, std::string triggerCh, std::string extra, std::string primary)
{
    int out = 0;
    std::string mkdir = "mkdir -p ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    out = system(mkdir.c_str());
    std::string folder = folder_name(run,subrun,voltage,threshold,triggerCh);

    std::string mv0 = "mv -n ~/Desktop/WaveDumpData/wave0.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    std::string mv1 = "mv -n ~/Desktop/WaveDumpData/wave1.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    std::string mv2 = "mv -n ~/Desktop/WaveDumpData/wave2.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    std::string mv3 = "mv -n ~/Desktop/WaveDumpData/wave3.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";

    std::string voltageS = changeVoltage(voltage);

    //QMessageBox::about(this,"",QString::fromStdString(folder0));
    mkdir = mkdir+folder;
    mv0 = mv0 + folder;
    mv1 = mv1 + folder;
    mv2 = mv2 + folder;
    mv3 = mv3 + folder;

    // checking if there is no space in the extra
    char extraChar[extra.size()+1];
    strcpy(extraChar,extra.c_str());
    char c;
    int aux = 0;
    bool noSpace = true;
    while(extraChar[aux]){
         c = extraChar[aux];
         if(isspace(c)){
             noSpace = false;
         }
         aux++;
    }

    // checked

    mv0 = mv0 + std::to_string(subrun) + "_wave0_" + voltageS + "_" + std::to_string(static_cast<int>(threshold)) + "ADC_" + triggerCh;
    mv1 = mv1 + std::to_string(subrun) + "_wave1_" + voltageS + "_" + std::to_string(static_cast<int>(threshold)) + "ADC_" + triggerCh;
    mv2 = mv2 + std::to_string(subrun) + "_wave2_" + voltageS + "_" + std::to_string(static_cast<int>(threshold)) + "ADC_" + triggerCh;
    mv3 = mv3 + std::to_string(subrun) + "_wave3_" + voltageS + "_" + std::to_string(static_cast<int>(threshold)) + "ADC_" + triggerCh;
    if(extra!=""){
        if(noSpace){
            mv0 = mv0 + "_" + extra;
            mv1 = mv1 + "_" + extra;
            mv2 = mv2 + "_" + extra;
            mv3 = mv3 + "_" + extra;
        }
        else{
            QMessageBox::about(this,"","Warning: extra has space");
            return false;
        }
    }
    mv0 = mv0+".dat";
    mv1 = mv1+".dat";
    mv2 = mv2+".dat";
    mv3 = mv3+".dat";

    out = system(mkdir.c_str());

    //QMessageBox::about(this,"",QString::fromStdString(mv0));
    out = system(mv0.c_str());
    out = system(mv1.c_str());
    out = system(mv2.c_str());
    out = system(mv3.c_str());


    return true;


}

std::string MainWindow::changeVoltage(double voltage){
    std::string post;
    std::string pre;
    pre = std::to_string(static_cast<int>(voltage));
    post = std::to_string(static_cast<int>(trunc(voltage*10))%10) + std::to_string(static_cast<int>(trunc(voltage*1000))%100 / 10);
    std::string voltageS = pre+"V"+post;
    //QMessageBox::about(this,"",QString::fromStdString(voltageS));

    return voltageS;
}


void MainWindow::on_button_movefile_2_clicked()
{
    QMessageBox::about(this,"","Calibration file moved. Check it!");

    ui->calibration_check->setChecked(true);
    // Take info from the data tab, so it is possible to go to the right folder
    int runNo = std::stoi(ui->run->text().toStdString());
    int subRunNo = std::stoi(ui->subrun->text().toStdString());
    double voltage = std::stod(ui->voltage->text().toStdString());
    double threshold = std::stod(ui->threshold->text().toStdString());
    std::string triggerCh = ui->trigger_channel->text().toStdString();
    std::string extra = ui->extra->text().toStdString();
    std::string primary = ui->primary_name->text().toStdString();

    move_calibration_file(runNo,subRunNo,voltage,threshold,triggerCh,extra,primary);

}

bool MainWindow::move_calibration_file(int run, int subrun, double voltage, double threshold, std::string triggerCh, std::string extra, std::string primary)
{
    int out = 0;
    std::string mkdir = "mkdir -p ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    out = system(mkdir.c_str());
    std::string folder = folder_name(run,subrun,voltage,threshold,triggerCh);
    folder = folder + "Calibration/";

    std::string mv0 = "mv -n ~/Desktop/WaveDumpData/wave0.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    std::string mv1 = "mv -n ~/Desktop/WaveDumpData/wave1.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    std::string mv2 = "mv -n ~/Desktop/WaveDumpData/wave2.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";
    std::string mv3 = "mv -n ~/Desktop/WaveDumpData/wave3.dat ~/Documents/ADC_data/x_arapuca_argon_test/" + primary + "/";

    std::string voltageS = changeVoltage(voltage);

    double led_voltage = std::stod(ui->led_voltage->text().toStdString());
    std::string led_voltageS = changeVoltage(led_voltage);
    int width = std::stoi(ui->led_width->text().toStdString());

    //QMessageBox::about(this,"",QString::fromStdString(folder0));
    mkdir = mkdir+folder;
    mv0 = mv0 + folder;
    mv1 = mv1 + folder;
    mv2 = mv2 + folder;
    mv3 = mv3 + folder;

    mv0 = mv0 + "wave0_" + voltageS + "_" + led_voltageS + "_" + std::to_string(width) + "ns";
    mv1 = mv1 + "wave1_" + voltageS + "_" + led_voltageS + "_" + std::to_string(width) + "ns";
    mv2 = mv2 + "wave2_" + voltageS + "_" + led_voltageS + "_" + std::to_string(width) + "ns";
    mv3 = mv3 + "wave3_" + voltageS + "_" + led_voltageS + "_" + std::to_string(width) + "ns";

    mv0 = mv0+".dat";
    mv1 = mv1+".dat";
    mv2 = mv2+".dat";
    mv3 = mv3+".dat";

    out = system(mkdir.c_str());

    //QMessageBox::about(this,"",QString::fromStdString(mv0));
    out = system(mv0.c_str());
    out = system(mv1.c_str());
    out = system(mv2.c_str());
    out = system(mv3.c_str());

    return true;
}


void MainWindow::on_lock_folder_clicked(bool checked)
{
    if(checked){
        ui->primary_name->setDisabled(true);
    }
    else{
        ui->primary_name->setEnabled(true);
    }
}
