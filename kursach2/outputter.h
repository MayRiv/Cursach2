#ifndef OUTPUTTER_H
#define OUTPUTTER_H

#include <QDialog>
#include <QVector>
#include <QString>
#include <QTextStream>

#include <qwt_plot.h>
#include <qwt3d_surfaceplot.h>

#include "triplet.h"
namespace Ui {
class Outputter;
}

class Outputter : public QDialog
{
    Q_OBJECT
    
public:
    explicit Outputter(QWidget *parent = 0);
    ~Outputter();
    void addGraph(std::pair<QVector<double>,QVector<double> > _graph, int tabNumber,QColor color);
    void addGraph(QVector<double> x,QVector<double>  y,int tabNumber,QColor color);
    void addGraph(double* x, double* y, int numberOfElements, int tabNumber, QColor color);
    void addGraph(QVector<double> x,QVector<double> y, QVector<double> z,int tabNumber);
    void view();
    void printf(QString string);
    QTextStream stream;
private:
    template<class T>
    bool checkPlotDimension(int tabNumber,QVector<std::pair<T,int> > plots);
    int numberOfTabs;
    Ui::Outputter *ui;
    QVector<std::pair< std::pair <QVector <double> ,QVector <double> >,std::pair<int,QColor> > > graphs;
    QVector<std::pair<Triplet<QVector<double> >,int> > graphs3d;

    QVector<std::pair<QwtPlot*,int> > plots;
    QVector<std::pair<Qwt3D::SurfacePlot*,int> > plots3d;
    QString outputString;
};
#endif // OUTPUTTER_H
