#include "WdgtKinectSuperPoints.h"

#include <QtGui>
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    WdgtKinectSuperPoints w;
    w.show();
    return a.exec();
}
