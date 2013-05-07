
#include "calculator.h"
#include "ui_calculator.h"
#include <math.h>
#include <iostream>
Calculator::Calculator(QWidget *parent,Outputter* out, int _nX, int _nT) :
    QDialog(parent),
    ui(new Ui::Calculator)
{
    leftBoundary=0;
    rightBoundary=1;
    _out=out;
    ui->setupUi(this);
    double h=(rightBoundary-leftBoundary)/(_nX-1);
    x.push_back((QVector<double>)0);
    x[0].push_back(leftBoundary);
    for (int i=1;i<_nX;i++)
        x[0].push_back(x[0].back()+h);
    a=3;
   // if (t/pow(h,2)>1.0/6) exit(123);

    nX=_nX;
    nT=_nT;
}

Calculator::~Calculator()
{
    delete ui;
}

void Calculator::calculate()
{
    double t=(rightBoundary-leftBoundary)/(nT-1);
    double h=(rightBoundary-leftBoundary)/(nX-1);
    u.push_back((QVector<double>)0);
    z.push_back((QVector<double>)0);

    foreach(double node,x[0])
    {
        u[0].push_back(getAccurateValue(node,0));   //Beginning
        z[0].push_back(getAccurateValue(node,0));  //conditions
    }
    double time=leftBoundary;
    y.push_back(time);

    while(time<rightBoundary-t)
    {
        time+=t;
        y.push_back(time);
        z.push_back((QVector<double>)0);
        for (int i=0;i<x[0].size();i++)
        {
            z.back().push_back(getAccurateValue(x[0][i],time));
        }

        QVector<double> uTH(x.back().size());
        QVector<double> uSupport;
        QVector<double> oldX=x.back();
        QVector<double> erors;
        QVector<double> doubleX;
        QVector<double> uTHdiv2;
        QVector<double> uTdiv2H;
        QVector<double> uClarify;
        doubleX=getDoubleX(oldX);

        uSupport=solveInterpolation(oldX,u.back(),doubleX);
        uTHdiv2=calculateNewton(uSupport,time,h/2,t);              //STEP HAVE BE ANOTHER

        uTdiv2H=calculateNewton(u.back(),time-t/2,h,t/2);
        uTdiv2H=calculateNewton(uTdiv2H,time,h,t/2);

        uTH=calculateNewton(u.back(),time,h,t);

        double eps=getEps(uTH,uTdiv2H,uTHdiv2);
        uClarify=clarifyU(uTH,uTdiv2H,uTHdiv2);
        u.push_back(uClarify);//uTH);
        double test=getMax(x.back());

        x.push_back(createNewWeb(oldX, erors));

        double argument=leftBoundary;
        foreach (double value, u.back())
        {
            std::cout <<"time is " << time <<" yApp=  " << value <<" yAcc=  " <<getAccurateValue(argument,time) <<" absPoh=  "<< value-getAccurateValue(argument,time)  << " otnPoh=  "<<(value-getAccurateValue(argument,time))/value*100. <<"\n";
            argument+=h;
        }

        _out->stream << "\n";

    }

}

QVector<double> Calculator::calculateNewton(QVector<double> oldU,double time,double h, double t)
{
    QVector<double> result(oldU.size());
    result=oldU;

    result[0]=getAccurateValue(leftBoundary,time);
    result.back()=getAccurateValue(rightBoundary,time);

    QVector<double> du(result.size());
    QVector<double> A(du.size()*du.size());
    QVector<double> B(result.size());

    for (int i=0;i<du.size();i++)
        du[i]=5;


    while (fabs(du[3])>0.005)
    {
        B[0]=0;
        B[B.size()-1]=0;
        A=fillYacoby(result,h,t);
        for (int i=1;i<B.size()-1;i++)
            B[i]=-fi(oldU,result[i],result[i+1],result[i-1],i,h,t);
        du=solveGauss(A,B);
        for (int i=0;i<du.size();i++)
            result[i]=result[i]+du[i];

    }
    return result;
}

QVector<double> Calculator::fillYacoby(QVector<double> us,double h, double t)
{
    QVector<double> A(us.size()*us.size());
    A[0]=1;
    for (int i=1;i<us.size()-1;i++)
    {
        A[i*us.size()+i]=dfdui(us[i],us[i+1],us[i-1],h,t);
        A[i*us.size()+i+1]=dfduiplus1(us[i],us[i+1],us[i-1],h,t);
        A[i*us.size()+i-1]=dfduiminus1(us[i],us[i+1],us[i-1],h,t);
    }
    A[A.size()-1]=1;
    return A;
}

QVector<double> Calculator::createNewWeb(QVector<double> oldX, QVector<double> erors)
{
    return oldX;   //temporary
}
QVector<double> Calculator::getDoubleX(QVector<double> oldX)
{
    QVector<double> newX(oldX.size()*2);
    double h=(rightBoundary-leftBoundary)/(newX.size()-1);
    newX[0]=leftBoundary;
    for (int i=1;i<newX.size();i++)
        newX[i]=newX[i-1]+h;
    return newX;
}
double  Calculator::getEps(QVector<double> uTH, QVector<double> uTdiv2H, QVector<double> uTHdiv2)
{
    QVector<double> eps(uTH.size());
    for (int i=0;i<eps.size();i++)
    {
        eps[i]=2 * uTdiv2H[i]+4.0/3 * uTHdiv2[i*2]-10.0/3 * uTH[i];
    }
    double maxEps=getMax(eps);
    return maxEps;
}

QVector<double> Calculator::clarifyU(QVector<double> uTH, QVector<double> uTdiv2H, QVector<double> uTHdiv2)
{
    QVector<double> clU(uTH.size());
    for (int i=0;i<clU.size();i++)
    {
        clU[i]=2 * uTdiv2H[i]+4.0/3 * uTHdiv2[i*2]-7.0/3 * uTH[i];
    }
    return clU;
}

double Calculator::dfdui(double ui, double uiplus1, double uiminus1,double h, double t)
{

    double sigma=2*a*t/pow(2*h,2);                     //2 is for implicide method, without 2 is for Krank-Nikolson;
    double ksi=2*a*t/(2*h*h);
    return -1 + 2*ksi*ui*(-3*ui+uiplus1+uiminus1)+sigma*pow(uiplus1-uiminus1,2);
}

double Calculator::dfduiplus1(double ui,double uiplus1, double uiminus1,double h, double t)
{
    double sigma=2*a*t/pow(2*h,2);
    double ksi=2*a*t/(2*h*h);
    return ui*(ui*ksi+2*sigma*(uiplus1-uiminus1));
}

double Calculator::dfduiminus1(double ui,double uiplus1, double uiminus1,double h, double t)
{

    double sigma=2*a*t/pow(2*h,2);
    double ksi=2*a*t/(2*h*h);
    return ui*(ui*ksi+2*sigma*(uiminus1-uiplus1));
}
double Calculator::fi(QVector<double> oldU,double ui, double uiplus1, double uiminus1,int i,double h, double t)
{
    double sigma=2*a*t/pow(2.0*h,2);
    double ksi=2*a*t/(2.0*h*h);
    return oldU[i]-ui+sigma*ui*pow((uiplus1-uiminus1),2)+ksi*pow(ui,2)*(uiminus1-2*ui+uiplus1);
            //+sigma*oldU[i]*pow((oldU[i+1]-oldU[i-1]),2)+ksi*pow(oldU[i],2)*(oldU[i-1]-2*oldU[i]+oldU[i+1]);  This string is for Krank-Nikolson

}

double Calculator::getAccurateValue(double x, double y)
{
    double B;
    double A=B=3;
    return pow(pow((x-A),2)/(4*a*(B-y)),1.0/2);

}
double* Calculator::methodGauss02(const double* pA,	const double* pB,	int n )
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
QVector<double> Calculator::solveGauss(QVector<double> A, QVector<double> B)//Векторная обёртка над массивным гауссом.
{
   double* pA=new double [A.size()];
   double* pB=new double [B.size()];
   for (int i=0;i<A.size();i++) pA[i]=A[i];
   for (int i=0;i<B.size();i++) pB[i]=B[i];
   double* pX=new double [B.size()];
   pX=methodGauss02(pA,pB,B.size());
   QVector<double> result(B.size());
   for (int i=0;i<B.size();i++) result[i]=pX[i];
   delete pA;
   delete pB;
   delete pX;
   return result;

}

QVector<double> Calculator::solveInterpolation(QVector<double> xOld, QVector<double> yOld, QVector<double> xNew)
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
double Calculator::getMax(QVector<double> array)
{
    double max=array[0];
    foreach(double value, array)
    {
        if (fabs(max)<fabs(value)) max=value;
    }
    return max;
}
