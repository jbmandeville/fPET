#include "simwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication qt_app(argc, argv);
    SimWindow win;

    const QIcon *positronIcon = new QIcon(":/My-Icons/retirement.png");
    qt_app.setWindowIcon(*positronIcon);

    win.show();
    return qt_app.exec();
}
