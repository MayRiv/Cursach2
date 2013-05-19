/********************************************************************************
** Form generated from reading UI file 'outputter.ui'
**
** Created: Sat Apr 13 13:27:01 2013
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OUTPUTTER_H
#define UI_OUTPUTTER_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QPlainTextEdit>
#include <QtGui/QTabWidget>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Outputter
{
public:
    QGridLayout *gridLayout;
    QTabWidget *tabWidget;
    QWidget *tabText;
    QGridLayout *gridLayout_2;
    QPlainTextEdit *plainTextEdit;

    void setupUi(QDialog *Outputter)
    {
        if (Outputter->objectName().isEmpty())
            Outputter->setObjectName(QString::fromUtf8("Outputter"));
        Outputter->resize(548, 511);
        gridLayout = new QGridLayout(Outputter);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        tabWidget = new QTabWidget(Outputter);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setEnabled(true);
        tabText = new QWidget();
        tabText->setObjectName(QString::fromUtf8("tabText"));
        tabText->setEnabled(true);
        gridLayout_2 = new QGridLayout(tabText);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        plainTextEdit = new QPlainTextEdit(tabText);
        plainTextEdit->setObjectName(QString::fromUtf8("plainTextEdit"));
        plainTextEdit->setReadOnly(true);

        gridLayout_2->addWidget(plainTextEdit, 0, 0, 1, 1);

        tabWidget->addTab(tabText, QString());

        gridLayout->addWidget(tabWidget, 0, 0, 1, 1);


        retranslateUi(Outputter);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(Outputter);
    } // setupUi

    void retranslateUi(QDialog *Outputter)
    {
        Outputter->setWindowTitle(QApplication::translate("Outputter", "Outputter", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tabText), QApplication::translate("Outputter", "Text", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class Outputter: public Ui_Outputter {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OUTPUTTER_H
