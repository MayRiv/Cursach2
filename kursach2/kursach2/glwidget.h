#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    explicit GLWidget(QWidget *parent = 0);
    void setPoints(QVector<double> arguments, QVector<double> values);
private:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    double getMax();
    QVector<std::pair<double,double> > points;
};

#endif // GLWIDGET_H
