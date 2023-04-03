#include "simwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication qt_app(argc, argv);
    SimWindow win;

    const QIcon *positronIcon = new QIcon(":/My-Icons/positron.png");
    qt_app.setWindowIcon(*positronIcon);

    QProcess process;
//    process.start("/bin/ls");
    process.startDetached("/bin/ls");
//    process.waitForFinished(); // sets current thread to sleep and waits for pingProcess end
//    QString output(process.readAllStandardOutput());
//    qInfo() << output;
//    QString DataAsString = QString(output);
//    qInfo() << DataAsString;

    win.show();
    return qt_app.exec();
}
