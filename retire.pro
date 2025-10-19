QT += gui printsupport core
DEFINES += QT_NO_DEBUG_OUTPUT

LIBS += -L"../util/3rdparty/nifticlib-2.0.0/znzlib" -L"../3rdparty/zlib-1.2.5/lib" -lznz -lz
LIBS += -L"../util/3rdparty/fftw-3.3.8/lib" -lfftw3

INCLUDEPATH += ../util
INCLUDEPATH += ../util//3rdparty/qcustomplot

VPATH += ../util
VPATH += ../util/3rdparty
VPATH += ../util/3rdparty/qcustomplot

CONFIG+=sdk_no_version_check
#CONFIG += c++11 console
#CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        main.cpp \
        plot.cpp \
        qcustomplot.cpp \
        simwindow_GUI.cpp \
        simwindow_engine.cpp \
        utilio.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    simwindow.h \
    plot.h \
    qcustomplot.h \
    io.h

DISTFILES +=

RESOURCES += \
    myresources.qrc
