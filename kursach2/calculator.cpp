
#include "calculator.h"
#include "ui_calculator.h"
#include <math.h>
#include <iostream>
#include <QDebug>
using namespace std;
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
    edop=0.01;
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
    QVector<double> doubleX;
    QVector<double> uTHdiv2;
    QVector<double> uTdiv2H;
    QVector<double> uClarify;
    QVector<double> oldX;
    QVector<double> uTH;
    QVector<double> uSupport;

    /*QVector<double> xs;
    xs.push_back(0);
    xs.push_back(3);
    xs.push_back(5);                 FOR A TEST
    xs.push_back(20);
    xs.push_back(40);
    QVector<double> s=getSteps(xs);*/

    while(time<rightBoundary-t)
    {
        time+=t;
        y.push_back(time);
        z.push_back((QVector<double>)0);
        for (int i=0;i<x[0].size();i++)
        {
            z.back().push_back(getAccurateValue(x[0][i],time));
        }

        uTH.resize(x.back().size());

        oldX=x.back();
        doubleX=getDoubleX(oldX);

        uSupport=solveInterpolation(oldX,u.back(),doubleX);
        uTHdiv2=calculateNewton(uSupport,time,t,getSteps(doubleX));

        uTdiv2H=calculateNewton(u.back(),time-t/2,t/2,getSteps(oldX));
        uTdiv2H=calculateNewton(uTdiv2H,time,t/2,getSteps(oldX));

        uTH=calculateNewton(u.back(),time,t,getSteps(oldX));

        double eps=getEps(uTH,uTdiv2H,uTHdiv2);
        if (fabs(eps)>edop )
        {

            time-=t;
            t/=2;
            if (t<tMin) t=tMin;
            QVector<double> supportX=getDoubleX(x.back());
            u.back()=solveInterpolation(x.back(),u.back(),supportX);
            x.back()=supportX;
            _out->stream << "FAIL\n";
        }
        else
        {
          uClarify=clarifyU(uTH,uTdiv2H,uTHdiv2);
          double alpha;
          QVector<double> bettas=getCoeffs(uClarify, uTdiv2H,uTHdiv2,getSteps(x.back()),t,alpha);

          x.push_back(createNewWeb(oldX, bettas,getSteps(oldX)));
          QVector<double> uFitsNewWeb=solveInterpolation(oldX,uClarify,x.back());
          u.push_back(uFitsNewWeb);
          t*=alpha;

          double argument=leftBoundary;
          foreach (double value, u.back())
          {
             _out->stream /*std::cout*/ <<"time is " << time <<" yApp=  " << value <<" yAcc=  " <<getAccurateValue(argument,time) <<" absPoh=  "<< value-getAccurateValue(argument,time)  << " otnPoh=  "<<(value-getAccurateValue(argument,time))/value*100. <<"\n";
             argument+=getSteps(x.back())[0];//h;
          }

         }
         _out->stream << "\n";

    }
 _out->stream << y.size();
}

QVector<double> Calculator::calculateNewton(QVector<double> oldU, double time, double t, QVector<double> steps)
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
        A=fillYacoby(result,oldU,steps,t);
        for (int i=1;i<B.size()-1;i++)
            B[i]=-fi(oldU,result[i],result[i+1],result[i-1],i,steps[i-1],steps[i],t);     //something weird with steps
        du=solveGauss(A,B);
        for (int i=0;i<du.size();i++)
            result[i]=result[i]+du[i];

    }
    return result;
}

QVector<double> Calculator::fillYacoby(QVector<double> us,QVector<double> oldU,QVector<double> h, double t)
{
    QVector<double> A(us.size()*us.size());
    A[0]=1;
    for (int i=1;i<us.size()-1;i++)
    {
        A[i*us.size()+i]=dfdui(us[i],us[i+1],us[i-1],h[i-1],h[i],t,oldU,i);
        A[i*us.size()+i+1]=dfduiplus1(us[i],us[i+1],us[i-1],h[i-1],h[i],t,oldU,i);
        A[i*us.size()+i-1]=dfduiminus1(us[i],us[i+1],us[i-1],h[i-1],h[i],t,oldU,i);
    }
    A[A.size()-1]=1;
    return A;
}

QVector<double> Calculator::createNewWeb(QVector<double> oldX, QVector<double> bettas, QVector<double> steps)
{

    double betta=getMin(bettas);
    double h=steps[0];
    h*=betta;
    if (h>hMax) h=hMax;
    int number=(rightBoundary-leftBoundary)/h+1;
    h=(rightBoundary-leftBoundary)/(number-1);//floor;
    QVector<double> x(number);
    x[0]=leftBoundary;
    for (int i=1;i<number;i++)
        x[i]=x[i-1]+h;
    //if (h*(number-1)!=rightBoundary) qDebug() << "error in a function create new web";
    if (x.back()!=rightBoundary)
    {
        qDebug() << "error in a function create new web";
    }
    //===============================================
    QList<double> bettaBodnya;
    bettaBodnya.push_front(0); //for Bodnya's magic
    for (int i=0;i<bettas.size();i++)
        bettaBodnya.push_back(bettas[i]);
    QList<double> xBodnya;
    for (int i=0;i<oldX.size();i++)
        xBodnya.push_back(oldX[i]);

    //===============================================
    //QList<double> test = calc_new_hx_CS(xBodnya,bettaBodnya);
    //QVector<double> test=getTestWeb(oldX,bettas);
    return x;
}
QVector<double> Calculator::getDoubleX(QVector<double> oldX)
{
    QVector<double> stepsOld=getSteps(oldX);
    QVector<double> stepsNew;
    for (int i=0;i<stepsOld.size();i++)
    {
        stepsNew.push_back(stepsOld[i]/2);
        stepsNew.push_back(stepsOld[i]/2);
    }

    QVector<double> newX(oldX.size()*2-1);
    newX[0]=leftBoundary;
    for (int i=1;i<newX.size();i++)
    {
        newX[i]=newX[i-1]+stepsNew[i-1];
    }
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
QVector<double> Calculator::getCoeffs(QVector<double> uTH, QVector<double> uTdiv2H, QVector<double> uTHdiv2,QVector<double> steps, double t,double& alphaOut)
{
    QVector<double> hc(steps.size()-1);
    for (int i=0;i<hc.size();i++)
        hc[i]=(steps[i]+steps[i+1])/2.0;

    QVector<double> C1(uTH.size()-2),C2(uTH.size()-2);
    for (int i=1;i< uTH.size()-1;i++)
    {
        C1[i-1]=2.0/pow(t,2)*(uTdiv2H[i]-uTH[i]);
        C2[i-1]=4.0/(3.0*pow(hc[i-1],2)*t)*(uTHdiv2[i*2]-uTH[i]);
    }

    QVector<double> alphai(uTH.size()-2),bettai(uTH.size()-2);
    for (int i=1;i< uTH.size()-1;i++)
    {
        alphai[i-1]=sqrt(edop/(2*fabs(C1[i-1])*pow(t,2)));
    }
    double alpha=getMin(alphai);
    QVector<double> e1(uTH.size()-2);
    for (int i=0;i<e1.size();i++)
     {
        e1[i]=edop-fabs(C1[i])*pow(alpha*t,2);
        bettai[i]=1.0/hc[i]*sqrt(e1[i]/(alpha*t*fabs(C2[i])));
     }
    /*for (int i=0;i< bettai.size();i++)
        if (bettai[i]>4.0) bettai[i]=4.0;*/
    if (alpha>4.0) alpha=4.0;
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
    double y=uiplus1;
    double x=ui;
    double z=uiminus1;
    double h=hi;
    double p=hp1;
    double c=(p+h)/2;
    //return (fi(uOld,ui,uiplus1+eps,uiminus1,i,hi,hp1,t)-fi(uOld,ui,uiplus1,uiminus1,i,hi,hp1,t))/eps;
    return a*t*x*(p*(c*x+p*(x-z))+h*h*(y-x))/(c*c*p*p);
}

double Calculator::dfduiminus1(double ui, double uiplus1, double uiminus1, double hi, double hp1, double t,QVector<double> uOld,int i)
{

    double sigma=2*a*t/pow(2*hi,2);
    double ksi=2*a*t/(2*hi*hi);


    double eps=0.0001;
    double y=uiplus1;
    double x=ui;
    double z=uiminus1;
    double h=hi;
    double p=hp1;
    double c=(p+h)/2;

    //return (fi(uOld,ui,uiplus1,uiminus1+eps,i,hi,hp1,t)-fi(uOld,ui,uiplus1,uiminus1,i,hi,hp1,t))/eps;
    return a*t*x*(c*h*x+h*h*(x-y)+p*p*(z-x))/(c*c*h*h);
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
QVector<double> Calculator::getSteps(QVector<double> x)
{
    QVector<double> steps;
    for (int i=0;i<x.size()-1;i++)
        steps.push_back(x[i+1]-x[i]);
    return steps;
}



/*
QList<double> Calculator::calc_new_hx_CS(QList<double> oldx, QList<double> BetaI)
{
    QList<double> smth;
    int i,j,k1;
    bool quit,limit;
    QList<double> XmasOld = oldx;
    QList<double> XmasNew, HmasNew;
    double x0,x1,x2,y1,y2,prevstep,step,hr,rjm1,rj,muw,mum,h1,h2,kf;
    int Nold = XmasOld.size()-1;
    //{------форм-е начального графика допустимых шагов----}
    cout << "Hdp starting\n";
    QList<double>Hdp;
    Hdp.append(XmasOld[1]-XmasOld[0]);
    for (int i=1;i<BetaI.size();i++)
    {
        Hdp.append(BetaI[i] * ((XmasOld[i+1]-XmasOld[i-1])/2.0) );
    }
    Hdp.append(XmasOld[Nold]-XmasOld[Nold-1]);
    cout << "hdp cacled\n";
    cout<<"HDP: ";
    for (int i=0;i<Hdp.size();i++)
        cout << Hdp[i] << " ";
    cout<<"\n";

    //{------корр-ка графика допустимых шагов----}
    muw=1.4; mum=0.7;
    for (int i=1;i<=Nold-1;i++)
    {
        h1=Hdp[i-1]; h2=Hdp[i]; kf=h2/h1;
        if (h2>=h1)
        {
            if (kf>muw) Hdp[i]=h1*muw;
        }
        else
        {
            if (kf<mum)
            {
                 Hdp[i-1]=h2/mum;
                 for (int j=i-1;j>=1;j--)
                 {
                      h1=Hdp[j-1]; h2=Hdp[j]; kf=h2/h1;
                      if (h2>=h1)
                      {
                         if (kf>muw) Hdp[j]=h1*muw;
                      }
                      else
                      {
                         if (kf<mum) Hdp[j-1]=h2/mum;
                      }
                 } //{of j};
            } //{движения назад};
        }//{ убывающей ф-ии};
    }//{движения вперед};

    cout<<"HDP (corrected): ";
    for (int i=0;i<Hdp.size();i++)
        cout << Hdp[i] << " ";
    cout<<"\n";

    //{------построение новой сетки----}
    XmasNew.append(XmasOld[0]);
    i = 0; j = 1;
    quit = false;//{quit=1,когда сетка вся построена}
    do
    {//{цикл по j,т.е. по всем  новым узлам сетки,выход- по quit}
        x1=XmasNew[j-1]; y1=h_dop(XmasOld,Hdp,x1); x0=x1;
        limit=false;//{limit=1,если очередной шаг построен}
        do
        {//{цикл по i,т.е.по малым отрезкам [x1,x2], выход- по limit}
            x2=XmasOld[i+1]; y2=h_dop(XmasOld,Hdp,x2);
            if (y2>y1) y2=y1;
            if (y2<=(x2-x0))
            {
            //{-----------------корень найден------------------}
                 step=(x0+y1-(y2-y1)/(x2-x1)*x1)/(1-(y2-y1)/(x2-x1))-x0;
                 if (step<hMin) step=hMin;
                 if (step>hMax) step=hMax;
                 if (j>1)
                 {
                     prevstep=(XmasNew[j - 1] - XmasNew[j - 2]);
                     hr=step/prevstep;
                     if (hr>muw) step=prevstep*muw;
                     else if (hr<mum) step=prevstep*mum;
                 }
                 if (step<hMin) step=hMin;
                 if (step>hMax) step=hMax;
                 HmasNew.append(step);
                 XmasNew.append(leftBoundary+step); limit=true;
                 j=j+1;
            }
            else if (i==Nold-1)
            { //{----последний шаг--------}
                cout << "последний шаг\n";
                 step=y2; HmasNew.append(step); XmasNew.append(x0+step);
                 limit=true; quit=true;
            }
            else
            {  // {------следующий малый отрезок-----------}
                x1=x2; y1=y2; i=i+1;
                cout << "следующий малый отрезок. x1 = "<<x1<<"; y1 = "<<y1<<";i = "<<i<<";Nold = "<<Nold<<"\n";
            }
        } while (limit!=true);
    } while (quit!=true);

#define DEBUG_STEP
#ifdef DEBUG_STEP
    cout << "вычисления закончены\n";
#endif
    QList <double> newh = HmasNew;
    QList <double> newx = XmasNew;

    cout << "j = "<<j<<"; XmasNew size-1 = "<<XmasNew.size()-1<<"\n";
    cout << "Hmasnew size-1 = "<<HmasNew.size()-1<<"\n";
    double LimitX = XmasOld.last();
    cout << "old lastx = "<< XmasOld.last()<<"\n";

    newh.clear();
    for (int i=0;i<newx.size()-1;i++)
        newh.append(newx[i+1] - newx[i]);
    if (fabs(newx.last() - rightBoundary) > 0.00001)
    {
#ifdef DEBUG_STEP
        cout << "превышена черта xn. РАСЧЕТ ЛЯМБДА\n";
#endif
        double lambda = rightBoundary / newx.last();
#ifdef DEBUG_STEP
        cout << "ЛЯМБДА = " << lambda << '\n';
        cout << "замена шагов\n";
        cout << "старые: ";
        for (int k=0;k<newh.size();k++) cout << newh[k] << ' ';
        cout << '\n';
#endif
        for (int k=0;k<newh.size();k++)
            newh[k] *= lambda;
#ifdef DEBUG_STEP
        cout << "новые : ";
        for (int k=0;k<newh.size();k++) cout << newh[k] << ' ';
        cout << '\n';


        cout << "замена исков\n";
        cout << "старые: ";
        for (int k=0;k<newx.size();k++) cout << newx[k] << ' ';
        cout << '\n';
#endif
        for (int k=0;k<newx.size()-1;k++)
            newx[k+1] = newx[k] + newh[k];

#ifdef DEBUG_STEP
        cout << "новые : ";
        for (int k=0;k<newx.size();k++) cout << newx[k] << ' ';
        cout << '\n';
#endif
    }

    QList<double> *result = new QList<double>;
    *result = newh;
    return *result;
}
*/
double Calculator::h_dop(QList<double> XmasOld, QList<double> Hdp, double x)
{
    int i;
    double res;
    bool quit;

    cout << "x = " << x << "\n";

    int Nold = XmasOld.size()-1;
    quit=false;
    i=0;
    while ((!quit) && (i<Nold))
    {
       if ((x>=XmasOld[i]) && (x<=XmasOld[i+1]))
       {
          quit=true;
          res=Hdp[i]+(x-XmasOld[i])/(XmasOld[i+1]-XmasOld[i])*(Hdp[i+1]-Hdp[i]);
       }
       i=i+1;
    }
    if (res<hMin) res=hMin;
    if (res>hMax) res=hMax;
    return res;

}
double Calculator::getHDop(QVector<double> oldX, QVector<double> betta, double x)
{
    QVector<double> steps=getSteps(oldX);

    if (x<leftBoundary+oldX[1])  return betta[0]*(steps[0]+steps[1])/2;
    if (x>rightBoundary-oldX[oldX.size()-2]) return betta.back()*(steps[steps.size()-2]+steps.back())/2;
    int i=1;
    while (!(x>oldX[i]&&x<oldX[i+1]))
        i++;
    return betta[i]*(steps[i]+steps[i+1])/2  +  (x-oldX[i])/steps[i+1]*(betta[i+1] * (steps[i+1]+steps[i+2])/2 -betta[i]*(steps[i]+steps[i+1])/2 );
}
QVector<double> Calculator::getTestWeb(QVector<double> oldX, QVector<double> betta)
{
    double dx=0.000005;
    double hMin=getHDop(oldX,betta,0.0);
    QVector<double> hNew;
    double sum=0;
    double xNew=0;
    double x=0;

    int j=1;
    while(sum<rightBoundary)
    {
        double hDop=getHDop(oldX,betta,x);
        if (hDop<hMin) hMin=hDop;
        if (x-xNew>hMin)
        {
            hNew.push_back(x-xNew);
            sum+=x-xNew;
            xNew=x-dx;
            j++;
        }
        x+=dx;
    }
    return hNew;
}
