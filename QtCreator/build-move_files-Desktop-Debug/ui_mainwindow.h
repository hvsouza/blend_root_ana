/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QLineEdit *primary_name;
    QTabWidget *tabWidget;
    QWidget *tab;
    QLabel *label_3;
    QLineEdit *run;
    QLineEdit *voltage;
    QPushButton *button_movefile;
    QLabel *label;
    QLineEdit *subrun;
    QPushButton *pushButton_2;
    QLineEdit *threshold;
    QLabel *label_6;
    QLabel *label_2;
    QLabel *label_5;
    QLabel *label_4;
    QLineEdit *extra;
    QLabel *label_7;
    QLabel *triggerCh;
    QLineEdit *trigger_channel;
    QWidget *tab_2;
    QLineEdit *led_width;
    QLabel *LED_voltage;
    QLineEdit *led_voltage;
    QLabel *width;
    QLabel *label_8;
    QLabel *label_9;
    QPushButton *button_movefile_2;
    QTabWidget *tabWidget_2;
    QWidget *tab_3;
    QWidget *tab_4;
    QLabel *Information;
    QLabel *label_10;
    QRadioButton *radioButton;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QMenuBar *menuBar;
    QMenu *menuLAr_Test;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(385, 372);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        MainWindow->setMinimumSize(QSize(385, 372));
        MainWindow->setMaximumSize(QSize(385, 372));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        primary_name = new QLineEdit(centralWidget);
        primary_name->setObjectName(QString::fromUtf8("primary_name"));
        primary_name->setEnabled(false);
        primary_name->setGeometry(QRect(161, 10, 211, 25));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(primary_name->sizePolicy().hasHeightForWidth());
        primary_name->setSizePolicy(sizePolicy1);
        primary_name->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        primary_name->setReadOnly(false);
        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setEnabled(true);
        tabWidget->setGeometry(QRect(10, 50, 361, 251));
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        label_3 = new QLabel(tab);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(24, 90, 67, 17));
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        run = new QLineEdit(tab);
        run->setObjectName(QString::fromUtf8("run"));
        run->setGeometry(QRect(94, 30, 91, 21));
        run->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        run->setReadOnly(false);
        voltage = new QLineEdit(tab);
        voltage->setObjectName(QString::fromUtf8("voltage"));
        voltage->setGeometry(QRect(94, 90, 91, 21));
        voltage->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        button_movefile = new QPushButton(tab);
        button_movefile->setObjectName(QString::fromUtf8("button_movefile"));
        button_movefile->setGeometry(QRect(250, 50, 89, 25));
        label = new QLabel(tab);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(24, 30, 67, 17));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        subrun = new QLineEdit(tab);
        subrun->setObjectName(QString::fromUtf8("subrun"));
        subrun->setGeometry(QRect(94, 60, 91, 21));
        subrun->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        pushButton_2 = new QPushButton(tab);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));
        pushButton_2->setGeometry(QRect(250, 140, 89, 25));
        QPalette palette;
        QBrush brush(QColor(0, 0, 0, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::WindowText, brush);
        QBrush brush1(QColor(239, 41, 41, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Button, brush1);
        QBrush brush2(QColor(255, 51, 51, 255));
        brush2.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Light, brush2);
        QBrush brush3(QColor(229, 25, 25, 255));
        brush3.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Midlight, brush3);
        QBrush brush4(QColor(102, 0, 0, 255));
        brush4.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Dark, brush4);
        QBrush brush5(QColor(136, 0, 0, 255));
        brush5.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Mid, brush5);
        palette.setBrush(QPalette::Active, QPalette::Text, brush);
        QBrush brush6(QColor(255, 255, 255, 255));
        brush6.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::BrightText, brush6);
        palette.setBrush(QPalette::Active, QPalette::ButtonText, brush);
        palette.setBrush(QPalette::Active, QPalette::Base, brush1);
        palette.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette.setBrush(QPalette::Active, QPalette::Shadow, brush);
        QBrush brush7(QColor(229, 127, 127, 255));
        brush7.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::AlternateBase, brush7);
        QBrush brush8(QColor(255, 255, 220, 255));
        brush8.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::ToolTipBase, brush8);
        palette.setBrush(QPalette::Active, QPalette::ToolTipText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::WindowText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Button, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Light, brush2);
        palette.setBrush(QPalette::Inactive, QPalette::Midlight, brush3);
        palette.setBrush(QPalette::Inactive, QPalette::Dark, brush4);
        palette.setBrush(QPalette::Inactive, QPalette::Mid, brush5);
        palette.setBrush(QPalette::Inactive, QPalette::Text, brush);
        palette.setBrush(QPalette::Inactive, QPalette::BrightText, brush6);
        palette.setBrush(QPalette::Inactive, QPalette::ButtonText, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Base, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Shadow, brush);
        palette.setBrush(QPalette::Inactive, QPalette::AlternateBase, brush7);
        palette.setBrush(QPalette::Inactive, QPalette::ToolTipBase, brush8);
        palette.setBrush(QPalette::Inactive, QPalette::ToolTipText, brush);
        palette.setBrush(QPalette::Disabled, QPalette::WindowText, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::Button, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Light, brush2);
        palette.setBrush(QPalette::Disabled, QPalette::Midlight, brush3);
        palette.setBrush(QPalette::Disabled, QPalette::Dark, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::Mid, brush5);
        palette.setBrush(QPalette::Disabled, QPalette::Text, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::BrightText, brush6);
        palette.setBrush(QPalette::Disabled, QPalette::ButtonText, brush4);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Shadow, brush);
        QBrush brush9(QColor(204, 0, 0, 255));
        brush9.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Disabled, QPalette::AlternateBase, brush9);
        palette.setBrush(QPalette::Disabled, QPalette::ToolTipBase, brush8);
        palette.setBrush(QPalette::Disabled, QPalette::ToolTipText, brush);
        pushButton_2->setPalette(palette);
        pushButton_2->setAutoFillBackground(false);
        pushButton_2->setStyleSheet(QString::fromUtf8("background-color: rgb(239, 41, 41);\n"
""));
        pushButton_2->setAutoDefault(false);
        pushButton_2->setDefault(false);
        pushButton_2->setFlat(false);
        threshold = new QLineEdit(tab);
        threshold->setObjectName(QString::fromUtf8("threshold"));
        threshold->setGeometry(QRect(94, 120, 91, 21));
        threshold->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_6 = new QLabel(tab);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(194, 120, 41, 16));
        label_2 = new QLabel(tab);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(24, 60, 67, 17));
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_5 = new QLabel(tab);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(194, 90, 16, 16));
        label_4 = new QLabel(tab);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(10, 120, 81, 20));
        label_4->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        extra = new QLineEdit(tab);
        extra->setObjectName(QString::fromUtf8("extra"));
        extra->setGeometry(QRect(94, 180, 91, 21));
        extra->setStyleSheet(QString::fromUtf8("background-color: rgb(211, 215, 207);"));
        extra->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_7 = new QLabel(tab);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(10, 180, 81, 20));
        label_7->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        triggerCh = new QLabel(tab);
        triggerCh->setObjectName(QString::fromUtf8("triggerCh"));
        triggerCh->setGeometry(QRect(10, 150, 81, 20));
        triggerCh->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        trigger_channel = new QLineEdit(tab);
        trigger_channel->setObjectName(QString::fromUtf8("trigger_channel"));
        trigger_channel->setGeometry(QRect(94, 150, 91, 21));
        trigger_channel->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 255, 255);"));
        trigger_channel->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        tabWidget->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        led_width = new QLineEdit(tab_2);
        led_width->setObjectName(QString::fromUtf8("led_width"));
        led_width->setGeometry(QRect(230, 50, 91, 21));
        led_width->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        LED_voltage = new QLabel(tab_2);
        LED_voltage->setObjectName(QString::fromUtf8("LED_voltage"));
        LED_voltage->setGeometry(QRect(136, 20, 91, 20));
        LED_voltage->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        led_voltage = new QLineEdit(tab_2);
        led_voltage->setObjectName(QString::fromUtf8("led_voltage"));
        led_voltage->setGeometry(QRect(230, 20, 91, 21));
        led_voltage->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        led_voltage->setReadOnly(false);
        width = new QLabel(tab_2);
        width->setObjectName(QString::fromUtf8("width"));
        width->setGeometry(QRect(136, 50, 91, 20));
        width->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_8 = new QLabel(tab_2);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(330, 20, 16, 16));
        label_9 = new QLabel(tab_2);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setGeometry(QRect(330, 50, 16, 16));
        button_movefile_2 = new QPushButton(tab_2);
        button_movefile_2->setObjectName(QString::fromUtf8("button_movefile_2"));
        button_movefile_2->setGeometry(QRect(180, 120, 161, 25));
        tabWidget_2 = new QTabWidget(tab_2);
        tabWidget_2->setObjectName(QString::fromUtf8("tabWidget_2"));
        tabWidget_2->setGeometry(QRect(30, 250, 127, 80));
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        tabWidget_2->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QString::fromUtf8("tab_4"));
        tabWidget_2->addTab(tab_4, QString());
        Information = new QLabel(tab_2);
        Information->setObjectName(QString::fromUtf8("Information"));
        Information->setGeometry(QRect(10, 20, 111, 151));
        Information->setAutoFillBackground(false);
        Information->setTextFormat(Qt::AutoText);
        Information->setScaledContents(false);
        Information->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        Information->setWordWrap(true);
        tabWidget->addTab(tab_2, QString());
        label_10 = new QLabel(centralWidget);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setGeometry(QRect(60, 10, 91, 20));
        label_10->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        radioButton = new QRadioButton(centralWidget);
        radioButton->setObjectName(QString::fromUtf8("radioButton"));
        radioButton->setGeometry(QRect(310, 40, 61, 23));
        radioButton->setChecked(true);
        MainWindow->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 385, 22));
        menuLAr_Test = new QMenu(menuBar);
        menuLAr_Test->setObjectName(QString::fromUtf8("menuLAr_Test"));
        MainWindow->setMenuBar(menuBar);

        menuBar->addAction(menuLAr_Test->menuAction());

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        primary_name->setText(QApplication::translate("MainWindow", "new_data", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "Voltage:", 0, QApplication::UnicodeUTF8));
        run->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        voltage->setText(QApplication::translate("MainWindow", "42.3", 0, QApplication::UnicodeUTF8));
        button_movefile->setText(QApplication::translate("MainWindow", "Move files", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Run:", 0, QApplication::UnicodeUTF8));
        subrun->setText(QApplication::translate("MainWindow", "0", 0, QApplication::UnicodeUTF8));
        pushButton_2->setText(QApplication::translate("MainWindow", "Finish run", 0, QApplication::UnicodeUTF8));
        threshold->setText(QApplication::translate("MainWindow", "20", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("MainWindow", "ADC", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "subrun:", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "V", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Threshold:", 0, QApplication::UnicodeUTF8));
        extra->setText(QString());
        label_7->setText(QApplication::translate("MainWindow", "Extra info:", 0, QApplication::UnicodeUTF8));
        triggerCh->setText(QApplication::translate("MainWindow", "Trigger Ch:", 0, QApplication::UnicodeUTF8));
        trigger_channel->setText(QApplication::translate("MainWindow", "Ch0", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("MainWindow", "Acquisition", 0, QApplication::UnicodeUTF8));
        led_width->setText(QApplication::translate("MainWindow", "100", 0, QApplication::UnicodeUTF8));
        LED_voltage->setText(QApplication::translate("MainWindow", "LED Voltage:", 0, QApplication::UnicodeUTF8));
        led_voltage->setText(QApplication::translate("MainWindow", "3", 0, QApplication::UnicodeUTF8));
        width->setText(QApplication::translate("MainWindow", "Width:", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("MainWindow", "V", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("MainWindow", "ns", 0, QApplication::UnicodeUTF8));
        button_movefile_2->setText(QApplication::translate("MainWindow", "Move Calibration files", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_3), QApplication::translate("MainWindow", "Tab 1", 0, QApplication::UnicodeUTF8));
        tabWidget_2->setTabText(tabWidget_2->indexOf(tab_4), QApplication::translate("MainWindow", "Tab 2", 0, QApplication::UnicodeUTF8));
        Information->setText(QApplication::translate("MainWindow", "If you are running Calibration before data acquisition, make sure that the Acquisition tab have the correct informations", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("MainWindow", "Calibration", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("MainWindow", "Folder name:", 0, QApplication::UnicodeUTF8));
        radioButton->setText(QApplication::translate("MainWindow", "Lock", 0, QApplication::UnicodeUTF8));
        menuLAr_Test->setTitle(QApplication::translate("MainWindow", "LAr Test", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
