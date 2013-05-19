#include "qfunc3d.h"
QFunc3D::QFunc3D(Qwt3D::SurfacePlot *sp):Function(sp){}
QFunc3D::QFunc3D():Function(){}
QFunc3D::QFunc3D(Qwt3D::SurfacePlot *sp,
    double *x,double *y,double *z,
    int N,int M) :
  Function(sp)
{
    setData(x,y,z,N,M);
}

QFunc3D::~QFunc3D() {
  free((void *)xd);
}

void QFunc3D::setData(double* x, double*y, double* z, int N, int M)
{
    xd = (double *)malloc((N+M+M*N)*sizeof(double));
    yd = xd + N;
    zd = yd + M;
    for (int k=0; k < N; k++) xd[k] = x[k];
    for (int q=0; q < M; q++)
    {
      yd[q] = y[q];
      for (int k=0; k < N; k++) zd[q*N+k] = z[q*N+k];
    }
    Nx = N; My = M;
}
void QFunc3D::setData(QVector<double> xVec,QVector<double>yVec,QVector<double> zVec)
{
    int N=xVec.size();
    int M=yVec.size();
    xd = (double *)malloc((N+M+M*N)*sizeof(double));
    yd = xd + N;
    zd = yd + M;
    for (int k=0; k < N; k++) xd[k] = xVec[k];
    for (int q=0; q < M; q++)
    {
      yd[q] = yVec[q];
      for (int k=0; k < N; k++) zd[q*N+k] = zVec[q*N+k];
    }
    Nx = N; My = M;
}
void QFunc3D::attach(Qwt3D::SurfacePlot *sp)
{
  Qwt3D::Function::assign(sp);
}
double QFunc3D::findMax(double* array,int size)
{
    double max=fabs(array[0]);
    for (int i=1;i<size;i++)
        if (fabs(array[i]) > max) max=fabs(array[i]);
    return max;
}
double QFunc3D::findMin(double * array, int size)
{
    double min = array[0];
    for(int i=0; i < size; i++)
        if (array[i] < min) min = array[i];
    return min;
}

void QFunc3D::setDomain()
{
    double maxX=findMax(xd,Nx);
    double minX=findMin(xd,Nx);
    double maxY=findMax(yd,My);
    double minY=findMin(yd,My);
    QFunc3D::Function::setDomain(minX,maxX,minY,maxY);
}

double QFunc3D::operator()(double x,double y)
{
  double dx = (xd[Nx-1]-xd[0]) / (Nx-1);
  int n = floor((x-xd[0])/dx + 0.5);
  if (n < 0 || n > Nx-1) return 0;
  double dy = (yd[My-1]-yd[0]) / (My-1);
  int m = floor((y-yd[0])/dy + 0.5);
  if (m < 0 || m > My-1) return 0;
  return zd[m*Nx+n];
}
