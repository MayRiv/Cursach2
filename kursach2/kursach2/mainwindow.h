#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "calculator.h"
#include <qwt-qt4/qwt_plot_curve.h>
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0, Outputter *out = 0);
    ~MainWindow();
    void setFunction(QVector<double> x, QVector<double> y);

private:
    Ui::MainWindow *ui;
    Outputter* _out;
    QVector<QVector<double> > _x,_u,_z;
    QVector<double> _y,_eps,_timeStep;
    QwtPlotCurve function;
private slots:

    void on_calcButton_clicked();
    void on_spinBoxLayers_valueChanged(int arg1);
    void drawNewGraphic();
    void on_comboBoxGraphs_currentIndexChanged(const QString &arg1);
    void on_comboBoxGraphs_textChanged(const QString &arg1);
};

#endif // MAINWINDOW_H
