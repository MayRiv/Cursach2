#include <QtGui/QApplication>
#include <QVector>
#include <math.h>
#include "calculator.h"
#include <QDialog>
#include "outputter.h"
double* methodGauss02(const double* pA,	const double* pB,	int n );
QVector<double> solveGauss(QVector<double> A, QVector<double> B);
QVector<double> solveProgonka(QVector<double> A, QVector<double> B);
void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x);
QVector<double> solveInterpolation(QVector<double> xOld, QVector<double> yOld, QVector<double> xNew);
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Outputter out;
    QDialog dialog;
    Calculator calculator(&dialog,&out,60,200);
    calculator.calculate();
    calculator.show();
    double A[16]={
        11,12, 0, 0,
        21,22,23, 0,
        0 ,32,33,34,
        0 ,0 ,43,44
    };
    QVector<double> arr(16);
    for (int i=0;i<16;++i)
        arr[i]=A[i];
    QVector<double> B(4);
    for (int i=0;i<4;i++)
        B[i]=i*i;
    //QVector<double> x1=solveProgonka(arr,B);
    //QVector<double> x2=solveGauss(arr,B);
    /*QVector<double> array(10),oldX(10),newX(20),newArray(20);
    array[0]=0;
    newX[0]=oldX[0]=0;
    double h=1.0/(10.0-1);
    int i;
    for (i=1;i<10;i++)
    {
        oldX[i]=oldX[i-1]+h;
        array[i]=sin(oldX[i]);
    }
    h/=2;
    for (int j=1;j<20;j++)
        newX[j]=newX[j-1]+h;
    newArray=solveInterpolation(oldX,array,newX);
    */
    out.view();
    out.show();
    return app.exec();
}
QVector<double> solveInterpolation(QVector<double> xOld, QVector<double> yOld, QVector<double> xNew)
        {
            const int sizeXNew = xNew.size();
            const int sizeXOld = xOld.size();
            QVector<double> yNew(sizeXNew);
            double Fi = 0;
            double p1 = 1;
            double p2 = 1;
            for (int k = 0; k < sizeXNew; k++)
            {
                Fi = 0;
                for (int i = 0; i < sizeXOld; i++)
                {
                    p1 = 1;
                    p2 = 1;
                    for (int j = 0; j < sizeXOld; j++)
                    {
                        if (i != j)
                        {
                            p1 = p1 * (xNew[k] - xOld[j]);
                            p2 = p2 * (xOld[i] - xOld[j]);
                        }
                    }
                    Fi = Fi + yOld[i] * p1 / p2;
                }
                yNew[k] = Fi;
            }
            return yNew;
        }
QVector<double> solveProgonka(QVector<double> A, QVector<double> B)
{
    double* above=new double [B.size()];
    double* diagonal=new double [B.size()-1];
    double* below=new double [B.size()-1];
    diagonal[0]=A[0];
    below[0]=A[B.size()];
    for (int i=1;i<B.size()-1;++i)
    {
        diagonal[i]=A[B.size()*i+i];
        below[i]=A[B.size()*(i+1)+i];
        above[i-1]=A[B.size()*(i-1)+i];
    }
    diagonal[B.size()-1]=A[B.size()*(B.size()-1)+(B.size()-1)];
    above[B.size()-2]=A[B.size()*(B.size()-2)+(B.size()-1)];
    QVector<double> first(B.size());
    QVector<double> second(B.size()-1);
    QVector<double> third(B.size()-1);
    for (int i=0;i<B.size()-1;++i)
    {
        first[i]=diagonal[i];
        second[i]=above[i];
        third[i]=below[i];
    }
    first[B.size()-1]=diagonal[B.size()-1];
    second[B.size()-2]=above[B.size()-2];

    double* f=new double [B.size()];
    for (int i=0;i<B.size();++i)
        f[i]=B[i];

    double* x=new double [B.size()];
    solveMatrix(B.size(),above,diagonal,below,f,x);
    QVector<double> result(B.size());
    for (int i=0;i<B.size();++i)
        result[i]=x[i];
    return result;
}
void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
        double m;
        for (int i = 1; i < n; i++)
        {
                m = a[i]/c[i-1];
                c[i] = c[i] - m*b[i-1];
                f[i] = f[i] - m*f[i-1];
        }

        x[n-1] = f[n-1]/c[n-1];

        for (int i = n - 2; i >= 0; i--)
                x[i]=(f[i]-b[i]*x[i+1])/c[i];

}
double* methodGauss02(const double* pA,	const double* pB,	int n )
{
    static const char* methodGauss02 = "methodGauss02!";
    static const char* nullPointer = "Null-pointer!";
    static const char* noMemory = "No-memory!";
    static const char* nEquals0 = "n = 0!";

    double* X = NULL;

    if(!pA || !pB)
        printf("\n%s: %s\n", methodGauss02, nullPointer);
    else if(!n)
        printf("\n%s: %s\n", methodGauss02, nEquals0);
    else
    {
        double** A = NULL;	//A[0] -> A[0][0]
        double* B = NULL;	//
        int* ci = NULL;	//current row indexes
        int* cj = NULL; //current column indexes

        if(!(A = (double**)malloc(sizeof(double*)*n)) ||
             !(B = (double*)malloc(sizeof(double)*n)) ||
             !(ci = (int*)malloc(sizeof(int)*n)) ||
             !(cj = (int*)malloc(sizeof(int)*n)) ||
             !(A[0] = (double*)malloc(sizeof(double)*n*n))
             )
             printf("\n%s: %s\n", methodGauss02, noMemory);
        else
        {
            if(!(X = (double*)malloc(sizeof(double)*n)))
                printf("\n%s: %s\n", methodGauss02, noMemory);
            else
            {
                int i;	//row index
                int j;	//column index
                int k;	//row index
                int imax;	//row index of maximum element
                int jmax;	//colon index of maximum element
                int temp;
                double R;	//R = A[i+1][i] / A[i][i]
                double S;	//Si = A[i][i+1] * X[i+1] + A[i][i+2] * X[i+2] -...

                for(i = 1; i < n; i++)
                    A[i] = &A[0][i*n];

                for(i = 0; i < n; i++)
                {
                    ci[i] = i;
                    cj[i] = i;
                }

                //Copying pA to A and B to B
                for(i = 0; i < n*n; i++)
                    A[0][i] = pA[i];
                for(i = 0; i < n; i++)
                    B[i] = pB[i];

                //Making A triangle
                for(i = 0; i < n-1; i++)
                {
                    //Seeking for maximum element
                    imax = i;
                    jmax = i;
                    for(k = i; k < n; k++)
                        for(j = i; j < n; j++)
                            if(fabs(A[ci[k]][cj[j]]) > fabs(A[ci[imax]][cj[jmax]]))
                            {
                                imax = k;
                                jmax = j;
                            }

                    //if imax != i then exchanging A[i] and A[imax] rows
                    if(imax != i)
                    {
                        temp = ci[i];
                        ci[i] = ci[imax];
                        ci[imax] = temp;
                    }
                    if(jmax != i)
                    {
                        temp = cj[i];
                        cj[i] = cj[jmax];
                        cj[jmax] = temp;
                    }

                    //Making 0s in i column from i+1 row to n-1 row
                    for(k = i+1; k < n; k++)
                    {
                        R = A[ci[k]][cj[i]] / A[ci[i]][cj[i]];

                        // for(j = i+1!!!...)
                        for(j = i; j < n; j++)
                            A[ci[k]][cj[j]] = A[ci[k]][cj[j]] - R * A[ci[i]][cj[j]];
                        B[ci[k]] = B[ci[k]] - R * B[ci[i]];
                    }
                }

                //Calculating X
                B[ci[n-1]] = B[ci[n-1]] / A[ci[n-1]][cj[n-1]];
                for(k = n-2; k >= 0; k--)
                {
                    S = 0;
                    for(j = k+1; j < n; j++)
                        S = S + A[ci[k]][cj[j]] * B[ci[j]];
                    B[ci[k]] = (B[ci[k]] - S) / A[ci[k]][cj[k]];
                }

                for(i = 0; i < n; i++)
                    X[cj[i]] = B[ci[i]];
            }
            free(A[0]);
        }

        free(A);
        free(B);
        free(ci);
        free(cj);
    }
    QVector<double> XTestVector(n);
    for (int i=0;i<n;i++) XTestVector[i]=X[i];
    return X;
}
QVector<double> solveGauss(QVector<double> A, QVector<double> B)//Векторная обёртка над массивным гауссом.
{
   double* pA=new double [A.size()];
   double* pB=new double [B.size()];
   for (int i=0;i<A.size();i++) pA[i]=A[i];
   for (int i=0;i<B.size();i++) pB[i]=B[i];
   double* pX=new double [B.size()];
   pX=methodGauss02(pA,pB,B.size());
   QVector<double> result(B.size());
   for (int i=0;i<B.size();i++) result[i]=pX[i];
   delete [] pA;
   delete [] pB;
   delete [] pX;
   return result;

}
