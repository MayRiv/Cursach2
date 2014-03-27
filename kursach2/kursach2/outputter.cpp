#include "outputter.h"
#include "ui_outputter.h"
#include <QLineEdit>
#include <qwt_plot.h>
#include <QLayout>
#include <QPen>
#include <qwt_plot_curve.h>
#include <algorithm>
#include <qfunc3d.h>
Outputter::Outputter(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Outputter)
{
    numberOfTabs=0;
    ui->setupUi(this);
    stream.setString(&outputString);
}

Outputter::~Outputter()
{
    delete ui;
}

void Outputter::addGraph(QVector<double> x,QVector<double>  y,int tabNumber, QColor color)
{   
    if (tabNumber>numberOfTabs+1) tabNumber=numberOfTabs+1;
    if (tabNumber<1) tabNumber=1;
    if (checkPlotDimension(tabNumber,plots3d)) tabNumber=numberOfTabs+1;
    graphs.push_back(std::pair< std::pair <QVector <double> ,QVector <double> >,std::pair<int,QColor> >(std::pair<QVector<double>,QVector<double> > (x,y),std::pair<int,QColor> (tabNumber,color)));
    if (tabNumber>numberOfTabs)
    {
      QwtPlot *plot=new QwtPlot;
      plots.push_back(std::pair<QwtPlot*, int> (plot,tabNumber));
      ui->tabWidget->addTab(plots.back().first,"Graph");
      numberOfTabs++;
    }
}
void Outputter::addGraph(QVector<double> x,QVector<double> y, QVector<double> z,int tabNumber)
{
    //SHOULD CHEK SIZE OF ALL VECTORS
    Triplet<QVector<double> > triplet(x,y,z);
    if (tabNumber>numberOfTabs+1) tabNumber=numberOfTabs+1;
    if (tabNumber<1) tabNumber=1;
    if (checkPlotDimension(tabNumber,plots)) tabNumber=numberOfTabs+1;
    graphs3d.push_back(std::pair<Triplet<QVector<double> >,int> (triplet,tabNumber));
    if (tabNumber>numberOfTabs)
    {
      Qwt3D::SurfacePlot* plot=new Qwt3D::SurfacePlot;
      plot->setPlotStyle(Qwt3D::FILLEDMESH);
      plots3d.push_back(std::pair<Qwt3D::SurfacePlot*,int>(plot,tabNumber));
      ui->tabWidget->addTab(plots3d.back().first,"Graph3d");
      numberOfTabs++;
    }
}
void Outputter::addGraph(double* x, double* y, int numberOfElements, int tabNumber, QColor color)
{
    QVector<double> xMas(numberOfElements);

    QVector<double> yMas(numberOfElements);

    for (int i=0;i<numberOfElements;i++)
    {
        xMas[i]=x[i];
        yMas[i]=y[i];
    }
    addGraph(xMas,yMas,tabNumber,color);
}
void Outputter::addGraph(std::pair<QVector<double>,QVector<double> > _graph,int tabNumber,QColor color)
{
    QVector<double> x=_graph.first;
    QVector<double> y=_graph.second;
    addGraph(x,y,tabNumber,color);
}

void Outputter::view()
{
    QVector<QwtPlotCurve*> functions(graphs.size());
    ui->plainTextEdit->clear();
    printf(outputString);
    for (int i=0;i<graphs.size();i++)
    {
        functions[i]=new QwtPlotCurve;
        functions[i]->setRenderHint(QwtPlotItem::RenderAntialiased);
        QPen pen(graphs[i].second.second);

        functions[i]->setPen(pen);

        functions[i]->setData(graphs[i].first.first,graphs[i].first.second);
        functions[i]->attach(plots[graphs[i].second.first-1].first);
    }

    for (int i=0;i<plots.size();i++)
    {
        plots[i].first->replot();
    }
    QVector<QFunc3D*> functions3d(graphs3d.size());
    for (int j=0;j<graphs3d.size();j++)
    {
        functions3d[j]=new QFunc3D;
        functions3d[j]->attach(plots3d[graphs3d[j].second-1].first);//here is a first problem - if we have mix of 3d graphs and 2d graphs, indexing to correct plot is wrong
        functions3d[j]->setData(graphs3d[j].first.first,graphs3d[j].first.second,graphs3d[j].first.third);
        functions3d[j]->setDomain();
        functions3d[j]->setMesh(graphs3d[j].first.first.size(),graphs3d[j].first.second.size());
        functions3d[j]->create();
    }

}

template<class T>
bool Outputter::checkPlotDimension(int tabNumber,QVector<std::pair<T,int> > plots)
{
    bool hasPlotDimension=false;
    for (int i=0;i<plots.size();i++)
    {
        if (tabNumber==plots[i].second) hasPlotDimension=true;
    }
    return hasPlotDimension;
}
void Outputter::printf(QString string)
{
    ui->plainTextEdit->setPlainText(string);
}
void Outputter::addGraph(QVector<double> x, QVector<double> y, QVector<QVector <double> > z, int tabNumber)
{
    QVector<double> A;
    foreach(QVector<double> vector,z)
        foreach(double value, vector)
            A.push_back(value);
    addGraph(x,y,A,tabNumber);
}
