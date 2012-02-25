#include "WdgtKinectSuperPoints.h"

#include <QtGui>
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    WdgtKinectSuperPoints w(argv[1]);
    w.show();
    return a.exec();
}
