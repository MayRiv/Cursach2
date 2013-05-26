#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>

class widget : public QWidget
{
    Q_OBJECT
public:
    explicit widget(QWidget *parent = 0);
    void setWeb(QVector<QVector<double> > x, QVector<double> t);
protected:
     void paintEvent(QPaintEvent *event);
    
private:
     QVector<std::pair <double,double > > function;
};

#endif // WIDGET_H
