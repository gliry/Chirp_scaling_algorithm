TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.c

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../fftw-3.3.10-build/ -llibfftw3.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../fftw-3.3.10-build/ -llibfftw3.dlld
else:unix: LIBS += -L$$PWD/../../../../fftw-3.3.10-build/ -llibfftw3.dll

INCLUDEPATH += $$PWD/../../../../fftw-3.3.10/api
DEPENDPATH += $$PWD/../../../../fftw-3.3.10/api
