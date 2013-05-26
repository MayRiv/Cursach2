#include "widget.h"
#include <QPainter>
widget::widget(QWidget *parent) :
    QWidget(parent)
{
    setFixedSize(1010, 710);
}

void widget::paintEvent(QPaintEvent *event)
{
    QPainter painter;
    painter.begin(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.drawPoint(50,40);

    for (int i=0;i<function.size();i++)
        painter.drawEllipse(1000*function[i].first,700-700*function[i].second,5,5);

    painter.end();
}
void widget::setWeb(QVector<QVector<double> > x, QVector<double> t)
{
    for (int i=0;i< t.size();i++)
    {
        foreach(double x, x[i])
            function.push_back(std::pair<double,double>(x,t[i]));
    }
}
