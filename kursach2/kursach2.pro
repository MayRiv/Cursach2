#-------------------------------------------------
#
# Project created by QtCreator 2013-04-22T14:45:18
#
#-------------------------------------------------

QT       += core gui

TARGET = Kursach2
TEMPLATE = app


SOURCES += main.cpp\
        calculator.cpp \
    qfunc3d.cpp \
    outputter.cpp

HEADERS  += calculator.h \
    qfunc3d.h \
    outputter.h \
    triplet.h


FORMS    += calculator.ui
unix {
INCLUDEPATH += /usr/include/qwt-qt4\
            /usr/include/qwtplot3d-qt4\
            /usr/include/GL
LIBS += -L/usr/lib -lqwt-qt4\
        -L/usr/lib -lqwtplot3d-qt4\
        -L/usr/lib/i386-linux-gnu -lGLU
}
