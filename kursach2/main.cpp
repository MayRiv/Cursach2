#include <QtGui/QApplication>
#include <QVector>
#include <math.h>
#include "calculator.h"
#include <QDialog>
#include "outputter.h"
QVector<double> solveInterpolation(QVector<double> xOld, QVector<double> yOld, QVector<double> xNew);
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    Outputter out;
    QDialog dialog;
    Calculator calculator(&dialog,&out,5,100);
    calculator.calculate();
    calculator.show();



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
