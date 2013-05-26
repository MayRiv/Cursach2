#include "glwidget.h"
#include <math.h>
GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent)
{
}
void GLWidget::initializeGL(){
    glClearColor(1,0,0,0);
}

void GLWidget::resizeGL(int w, int h){
    if (points.size()!=0)
        update();
}

void GLWidget::paintGL(){
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1,0,0);
    glBegin(GL_QUADS);
    glVertex3f(-1,-1,1);
    glVertex3f(1,-1,1);
    glVertex3f(1,1,1);
    glVertex3f(-1,1,1);
    glVertex3f(-1,-1,1);
    glEnd();
    glColor3f(0,1,1);
    double max;//=getMax()+0.1;
    glBegin(GL_LINE_STRIP);

    for (int i=0;i<points.size();i++)
    {
        glVertex3f(points[i].first/max,points[i].second/max,1);

    }

    glEnd();
    glBegin(GL_LINES);

    //axes
    glVertex2f(-0,-1);
    glVertex2f(-0,1);

    glVertex2f(-0,1);
    glVertex2f(-0.02,0.97);

    glVertex2f(-0,1);
    glVertex2f(0.02,0.97);


    glVertex2f(-1,0);
    glVertex2f(1,0);

    glVertex2f(1,0);
    glVertex2f(0.97,-0.03);
    glVertex2f(1,0);
    glVertex2f(0.97,0.03);
    glEnd();
}

void GLWidget::setPoints(QVector<double> arguments, QVector<double> values)
{
    if (arguments.size()!=values.size()) exit(40);
    for(int i=0;i<arguments.size();i++)
        points.push_back(std::pair<double,double>(arguments[i],values[i]));
}
double GLWidget::getMax(){
    double max=fabs(points[0].first);
    for (int i=1;i<points.size();i++)
        if (max<fabs(points[i].first)) max=points[i].first;
    for (int i=0;i<points.size();i++)
        if (max<fabs(points[i].second)) max=points[i].second;
    return max;
}
