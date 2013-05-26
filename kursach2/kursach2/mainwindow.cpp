#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "calculator.h"

#include <QPen>
MainWindow::MainWindow(QWidget *parent,Outputter* out) :
    QMainWindow(parent),_out(out),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setGeometry(0,0,1300,700);
    ui->setupUi(this);
    connect(ui->calcButton,SIGNAL(clicked()),this,SLOT(on_calcButton_clicked()));
    ui->comboBoxGraphs->addItem("Time step");
    ui->comboBoxGraphs->addItem("Real local error");
    ui->comboBoxGraphs->addItem("Difference between results");
    ui->implicit->setChecked(true);
    ui->spinBoxLayers->setRange(0,50);
    function.setRenderHint(QwtPlotItem::RenderAntialiased);
    function.attach(ui->qwtPlot);
    QPen pen(Qt::red);
    function.setPen(pen);

}

MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::setFunction(QVector<double> x, QVector<double> y)
{
    ui->widget->setPoints(x,y);
}

void MainWindow::on_calcButton_clicked()
{

    Calculator _c(this,_out,60,200);
    if (ui->implicit->isChecked()) _c.calculateImplicit();
    else    _c.calculateExplicit();
    _x=_c.xOut;//_c.getX();
    _y=_c.getTime();
    _u=_c.uOut;//_c.getU();
    _z=_c.getZ();
    _eps=_c.getEps();
    _timeStep=_c.timeStep;
    //ui->widget->setPoints(_y,_eps);
    //ui->widget->update();
    ui->widget->update();
    function.setData(_y,_eps);
    ui->qwtPlot->replot();
    ui->spinBoxLayers->setRange(0,_y.size()-1);
    drawNewGraphic();

}

void MainWindow::on_spinBoxLayers_valueChanged(int arg1)
{
    drawNewGraphic();
}
void MainWindow::drawNewGraphic()
{
    if (ui->comboBoxGraphs->currentText()=="Time step")
    {
        QVector<double> onetwothree(_timeStep.size());
        for(int i=0;i<_timeStep.size();i++)
            onetwothree[i]=i;
        function.setData(onetwothree,_timeStep);
    }
    if (ui->comboBoxGraphs->currentText()=="Real local error")
       {
        QVector<double> time;
        time.push_back(0);
        for(int i=0;i<_y.size();i++) time.push_back(time.back()+_y[i]);
        function.setData(time,_eps);
       }
    if (ui->comboBoxGraphs->currentText()=="Difference between results")
    {
        ui->spinBoxLayers->setEnabled(true);
        QVector<double> difference;
        int index=ui->spinBoxLayers->value();
        for(int i=0;i<_u[ui->spinBoxLayers->value()].size();i++)
            difference.push_back(fabs(_u[ui->spinBoxLayers->value()][i]-_z[ui->spinBoxLayers->value()][i]));
        function.setData(_x[ui->spinBoxLayers->value()],difference);
    }
    ui->qwtPlot->replot();
}




void MainWindow::on_comboBoxGraphs_currentIndexChanged(const QString &arg1)
{
    drawNewGraphic();
}

void MainWindow::on_comboBoxGraphs_textChanged(const QString &arg1)
{
    drawNewGraphic();
}
