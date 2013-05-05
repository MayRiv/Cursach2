#ifndef QFUNC3D_H
#define QFUNC3D_H

#include <qwt3d_function.h>

class QFunc3D : public Qwt3D::Function
{
public:
  QFunc3D(Qwt3D::SurfacePlot *,
    double *,double *,double *,
    int,int);
  QFunc3D(Qwt3D::SurfacePlot *sp);
  QFunc3D();
  ~QFunc3D();

  double operator ()(double,double);
  void attach(Qwt3D::SurfacePlot *sp);
  void setData(double* x, double*y, double* z, int N, int M);
  void setData(QVector<double> xVec,QVector<double>yVec,QVector<double> zVec);
  void setDomain();
private:
  double findMax(double*, int);
  double findMin(double*, int);
  double *xd,*yd,*zd;
  int Nx,My;
};

#endif // QFUNC3D_H
