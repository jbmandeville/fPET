QT += gui printsupport core
DEFINES += QT_NO_DEBUG_OUTPUT

LIBS += -L"../util/3rdparty/nifticlib-2.0.0/znzlib" -L"../3rdparty/zlib-1.2.5/lib" -lznz -lz

INCLUDEPATH += ../util
INCLUDEPATH += ../util//3rdparty/qcustomplot
INCLUDEPATH += ../FM
INCLUDEPATH += ../util/3rdparty/nifticlib-2.0.0/include
INCLUDEPATH += ../util/3rdparty/nifticlib-2.0.0/znzlib

VPATH += ../util
VPATH += ../FM
VPATH += ../util/3rdparty
VPATH += ../util/3rdparty/qcustomplot
VPATH += ../util/3rdparty/nifticlib-2.0.0/include

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
        simengine.cpp \
        simwindow.cpp \
        fmriglm.cpp \
        simanalyzer.cpp \
        generalglm.cpp \
        plot.cpp \
        qcustomplot.cpp \
        utilio.cpp \
        ImageIO.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    simengine.h \
    simwindow.h \
    fmriglm.h \
    simanalyzer.h \
    generalglm.h \
    plot.h \
    qcustomplot.h \
    ImageIO.h \
    io.h

DISTFILES +=

RESOURCES += \
    myresources.qrc
