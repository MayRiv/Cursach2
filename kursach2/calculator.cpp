
#include "calculator.h"
#include "ui_calculator.h"
#include <math.h>
#include <iostream>
#include <QDebug>
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
    a=5;
   // if (t/pow(h,2)>1.0/6) exit(123);
    edop=0.001;
    hMax=(rightBoundary-leftBoundary)/4;
    hMin=(rightBoundary-leftBoundary)/1000;
    tMin=(rightBoundary-leftBoundary)/1000;
    tMax=(rightBoundary-leftBoundary)/4;
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
       // QVector<double> betta;
        QVector<double> doubleX;
        QVector<double> uTHdiv2;
        QVector<double> uTdiv2H;
        QVector<double> uClarify;
        doubleX=getDoubleX(oldX);

        uSupport=solveInterpolation(oldX,u.back(),doubleX);
        uTHdiv2=calculateNewton(uSupport,time,h/2,t);

        uTdiv2H=calculateNewton(u.back(),time-t/2,h,t/2);
        uTdiv2H=calculateNewton(uTdiv2H,time,h,t/2);

        uTH=calculateNewton(u.back(),time,h,t);

        double eps=getEps(uTH,uTdiv2H,uTHdiv2);
        if (fabs(eps)>edop )
        {
            h/=2;
            if (h<hMin) h=hMin;
            time-=t;
            t/=2;
            if (t<tMin) t=tMin;
            QVector<double> supportX=getDoubleX(x.back());
            u.back()=solveInterpolation(x.back(),u.back(),supportX);
            x.back()=supportX;

        }
        else
        {
          uClarify=clarifyU(uTH,uTdiv2H,uTHdiv2);
          double alpha;
          QVector<double> bettas=getCoeffs(uClarify, uTdiv2H,uTHdiv2,h,t,alpha);

          x.push_back(createNewWeb(oldX, bettas,h));
          QVector<double> uFitsNewWeb=solveInterpolation(oldX,uClarify,x.back());
          //u.push_back(uClarify);
          u.push_back(uFitsNewWeb);


          t*=alpha;


          double argument=leftBoundary;
          foreach (double value, u.back())
          {
             _out->stream /*std::cout*/ <<"time is " << time <<" yApp=  " << value <<" yAcc=  " <<getAccurateValue(argument,time) <<" absPoh=  "<< value-getAccurateValue(argument,time)  << " otnPoh=  "<<(value-getAccurateValue(argument,time))/value*100. <<"\n";
             argument+=h;
          }

         }
         _out->stream << "\n";

    }
 _out->stream << y.size();
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


    while (fabs(getMax(du))>0.005)
    {
        B[0]=0;
        B[B.size()-1]=0;
        A=fillYacoby(result,oldU,h,t);
        for (int i=1;i<B.size()-1;i++)
            B[i]=-fi(oldU,result[i],result[i+1],result[i-1],i,h,h,t);
        du=solveGauss(A,B);
        for (int i=0;i<du.size();i++)
            result[i]=result[i]+du[i];

    }
    return result;
}

QVector<double> Calculator::fillYacoby(QVector<double> us,QVector<double> oldU,double h, double t)
{
    QVector<double> A(us.size()*us.size());
    A[0]=1;
    for (int i=1;i<us.size()-1;i++)
    {
        A[i*us.size()+i]=dfdui(us[i],us[i+1],us[i-1],h,h,t,oldU,i);
        A[i*us.size()+i+1]=dfduiplus1(us[i],us[i+1],us[i-1],h,h,t,oldU,i);
        A[i*us.size()+i-1]=dfduiminus1(us[i],us[i+1],us[i-1],h,h,t,oldU,i);
    }
    A[A.size()-1]=1;
    return A;
}

QVector<double> Calculator::createNewWeb(QVector<double> oldX, QVector<double> bettas,double& h)
{

    double betta=getMin(bettas);

    h*=betta;
    if (h>hMax) h=hMax;
    int number=(rightBoundary-leftBoundary)/h+1;
    h=(rightBoundary-leftBoundary)/(number-1);//floor;
    QVector<double> x(number);
    x[0]=leftBoundary;
    for (int i=1;i<number;i++)
        x[i]=x[i-1]+h;
    //x.last()=1;//kostul merzkiy                     //vse hernya tyt
    if (h*(number-1)!=rightBoundary) qDebug() << "error in a function create new web";
    if (x.back()!=rightBoundary)
    {
        qDebug() << "error in a function create new web";
       // exit(4);
    }
    //return oldX;
    return x;
}
QVector<double> Calculator::getDoubleX(QVector<double> oldX)
{
    QVector<double> newX(oldX.size()*2-1);
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
QVector<double> Calculator::getCoeffs(QVector<double> uTH, QVector<double> uTdiv2H, QVector<double> uTHdiv2,double h, double t,double& alphaOut)
{

    QVector<double> C1(uTH.size()),C2(uTH.size());
    for (int i=0;i< uTH.size();i++)
    {
        C1[i]=2.0/pow(t,2)*(uTdiv2H[i]-uTH[i]);
        C2[i]=4.0/(3.0*h*h*t)*(uTHdiv2[i*2]-uTH[i]);

    }

    QVector<double> alphai(uTH.size()),bettai(uTH.size());
    for (int i=0;i< uTH.size();i++)
    {
        alphai[i]=sqrt(edop/(2*fabs(C1[i])*pow(t,2)));
    }
    double alpha=getMin(alphai);
    QVector<double> e1(uTH.size());
    for (int i=0;i<e1.size();i++)
     {
        e1[i]=edop-fabs(C1[i])*pow(alpha*t,2);
        bettai[i]=1.0/h*sqrt(e1[i]/(alpha*t*fabs(C2[i])));
     }
    alphaOut=alpha;
    return bettai;
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

double Calculator::dfdui(double ui, double uiplus1, double uiminus1,double hi,double hp1, double t,QVector<double> uOld,int i)
{

    double sigma=2*a*t/pow(2*hi,2);                     //2 is for implicide method, without 2 is for Krank-Nikolson;
    double ksi=2*a*t/(2*hi*hi);
    //return -1 + 2*ksi*ui*(-3*ui+uiplus1+uiminus1)+sigma*pow(uiplus1-uiminus1,2);
    double y=uiplus1;
    double x=ui;
    double z=uiminus1;
    double h=hi;
    double p=hp1;
    double c=(p+h)/2;
    /*return a*t*(h*h*(3*x-y)+p*p*(z-3*x))*(h*h*(x-y)+p*p*(z-x))/(2*c*c*h*h*p*p)
           +
           a*t*x*(-3*c*x+h*y+p*z)/(c*h*p)
           -1;*/
    return 1/(2*c*c*h*h*p*p)*(-2*c*c*h*p*(6*a*t*x*x+h*p)+4*a*c*h*p*t*x*(h*y+p*z)+a*t*(pow(h,4)*(3*x*x-4*x*y+y*y)-2*h*h*p*p*(3*x*x-2*x*(y+z)+y*z)+pow(p,4)*(3*x*x-4*x*z+z*z) ) );
    double eps=0.000000001;
    //return (fi(uOld,ui+eps,uiplus1,uiminus1,i,hi,hp1,t)-fi(uOld,ui,uiplus1,uiminus1,i,hi,hp1,t))/eps;
}

double Calculator::dfduiplus1(double ui,double uiplus1, double uiminus1,double hi,double hp1, double t,QVector<double> uOld,int i)
{
    double sigma=2*a*t/pow(2*hi,2);
    double ksi=2*a*t/(2*hi*hi);
    //return ui*(ui*ksi+2*sigma*(uiplus1-uiminus1));
    double eps=0.0001;
    return (fi(uOld,ui,uiplus1+eps,uiminus1,i,hi,hp1,t)-fi(uOld,ui,uiplus1,uiminus1,i,hi,hp1,t))/eps;
}

double Calculator::dfduiminus1(double ui, double uiplus1, double uiminus1, double hi, double hp1, double t,QVector<double> uOld,int i)
{

    double sigma=2*a*t/pow(2*hi,2);
    double ksi=2*a*t/(2*hi*hi);
    //return ui*(ui*ksi+2*sigma*(uiminus1-uiplus1));
    double hc=(hp1+hi)/2;
    double x=ui;
    double z=uiplus1;
    double y=uiminus1;
    double h=hi;
    double p=hp1;
    double c=hc;
    double eps=0.0001;
    //return (fi(oldU,ui,uiplus1,uiminus1+eps,i,hi,hp1,t)-fi(oldU,ui,uiplus1,uiminus1,i,hi,hp1,t))/eps;
    //return a*t*x*(c*h*x+h*h*(x-z)+p*p*(y-z))/(c*c*h*h);

    return (fi(uOld,ui,uiplus1,uiminus1+eps,i,hi,hp1,t)-fi(uOld,ui,uiplus1,uiminus1,i,hi,hp1,t))/eps;

}

double Calculator::fi(QVector<double> oldU,double ui, double uiplus1, double uiminus1,int i,double hi,double hp1, double t)
{
    double hc=(hp1+hi)/2;
    double dudx=(uiplus1*pow(hi,2)+(pow(hp1,2)-pow(hi,2))*ui - uiminus1*pow(hp1,2))/(2*hc*hp1*hi);
    double d2udx2=(hp1 * uiminus1 -2 * hc * ui+ hi * uiplus1 ) / ( hc * hi * hp1 );
    double oldui=oldU[i];
    return oldui-ui+2*a*t*ui*pow(dudx,2)+a * t * pow(ui,2) * d2udx2;
    /*uiplus1=y;
    ui=x;
    uiminus1=z;
    hi=c;
    hp1=b;*/
}

double Calculator::getAccurateValue(double x, double y)
{
    double B;
    double A=B=3;
    //B=0.7;
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
   delete [] pA;
   delete [] pB;
   delete [] pX;
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
double Calculator::getMin(QVector<double> array)
{
    double min=array[0];
    foreach(double value, array)
    {
        if (fabs(min)>fabs(value)) min=value;
    }
    return min;
}
