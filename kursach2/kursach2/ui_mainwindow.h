/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Sun May 26 13:57:51 2013
**      by: Qt User Interface Compiler version 4.8.3
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qwt_plot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout;
    QSpacerItem *verticalSpacer;
    QRadioButton *implicit;
    QRadioButton *explicit_2;
    QComboBox *comboBoxGraphs;
    QSpinBox *spinBoxLayers;
    QPushButton *calcButton;
    QSpacerItem *verticalSpacer_2;
    QwtPlot *qwtPlot;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1140, 538);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout = new QGridLayout(centralwidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        implicit = new QRadioButton(centralwidget);
        implicit->setObjectName(QString::fromUtf8("implicit"));

        verticalLayout->addWidget(implicit);

        explicit_2 = new QRadioButton(centralwidget);
        explicit_2->setObjectName(QString::fromUtf8("explicit_2"));

        verticalLayout->addWidget(explicit_2);

        comboBoxGraphs = new QComboBox(centralwidget);
        comboBoxGraphs->setObjectName(QString::fromUtf8("comboBoxGraphs"));

        verticalLayout->addWidget(comboBoxGraphs);

        spinBoxLayers = new QSpinBox(centralwidget);
        spinBoxLayers->setObjectName(QString::fromUtf8("spinBoxLayers"));
        spinBoxLayers->setEnabled(false);

        verticalLayout->addWidget(spinBoxLayers);

        calcButton = new QPushButton(centralwidget);
        calcButton->setObjectName(QString::fromUtf8("calcButton"));

        verticalLayout->addWidget(calcButton);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer_2);


        gridLayout->addLayout(verticalLayout, 0, 0, 1, 1);

        qwtPlot = new QwtPlot(centralwidget);
        qwtPlot->setObjectName(QString::fromUtf8("qwtPlot"));

        gridLayout->addWidget(qwtPlot, 0, 1, 1, 1);

        MainWindow->setCentralWidget(centralwidget);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        implicit->setText(QApplication::translate("MainWindow", "Calc implicit", 0, QApplication::UnicodeUTF8));
        explicit_2->setText(QApplication::translate("MainWindow", "Calc explicit", 0, QApplication::UnicodeUTF8));
        calcButton->setText(QApplication::translate("MainWindow", "Calculate", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
