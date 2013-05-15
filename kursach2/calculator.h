#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <QDialog>
#include <QVector>
#include "outputter.h"
namespace Ui {
class Calculator;
}

class Calculator : public QDialog
{
    Q_OBJECT

public:
    explicit Calculator(QWidget *parent = 0,Outputter* out = 0,int nX=50, int nT=50);
    ~Calculator();
    void calculate();
private:
    double getAccurateValue(double x, double y);
    QVector<double> fillYacoby(QVector<double> us, QVector<double> oldU, QVector<double> h, double t);
    QVector<double> calculateNewton(QVector<double> oldU, double time, double t, QVector<double> steps);
    QVector<double> createNewWeb(QVector<double> oldX, QVector<double> bettas, QVector<double> steps);
    QVector<double> getDoubleX(QVector<double> oldX);
    double getEps (QVector<double> uTH, QVector<double> uTdiv2H, QVector<double> uTHdiv2);
    double getMax(QVector<double> array);
    double getMin(QVector<double> array);
    QVector<double> clarifyU(QVector<double> uTH, QVector<double> uTdiv2H, QVector<double> uTHdiv2);
    QVector<double> getCoeffs(QVector<double> uTH, QVector<double> uTdiv2H, QVector<double> uTHdiv2, QVector<double> steps,double t,double& alphaOut);
    QVector<double> getSteps(QVector<double> x);
    double* methodGauss02(const double* pA,	const double* pB,	int n );
    QVector<double> solveGauss(QVector<double> A, QVector<double> B);
    QVector<double> solveInterpolation(QVector<double> xOld, QVector<double> yOld, QVector<double> xNew);
    double dfdui(double ui, double uiplus1, double uiminus1, double hi, double hp1, double t, QVector<double> uOld, int i);
    double dfduiplus1(double ui, double uiplu1, double uiminus1, double hi, double hp1, double t, QVector<double> uOld, int i);
    double dfduiminus1(double ui, double uiplus1, double uiminus1, double hi, double hp1, double t, QVector<double> uOld, int i);
    double fi(QVector<double> oldU, double ui, double uiplus1, double uiminus1, int i, double hi, double hp1, double t);
    double getHDop(QVector<double> oldX, QVector<double> betta, double x);
    QVector<double> getTestWeb(QVector<double> oldX, QVector<double> betta);
    QList<double> calc_new_hx_CS(QList<double> oldx, QList<double> BetaI);
    double h_dop(QList<double> XmasOld, QList<double> Hdp, double x);
    QVector<double> y;//x=space, y=time
    QVector<QVector<double> > z,u,x;//z=accurate, u=approximate
    Ui::Calculator *ui;
    double a;
    int nX,nT;
    double leftBoundary,rightBoundary;
    double edop;
    double hMax,hMin;
    double tMin,tMax;
    Outputter* _out;
};

#endif // CALCULATOR_H
