#include "simwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication qt_app(argc, argv);
    SimWindow win;

    win.show();
    return qt_app.exec();
}
