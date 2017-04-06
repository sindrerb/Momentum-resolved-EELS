TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    kronigpenney.cpp \
    ../vec3-class/vec3.cpp \
    ../wavestate-class/wavestate.cpp

HEADERS += \
    kronigpenney.h \
    ../vec3-class/vec3.h \
    ../wavestate-class/wavestate.h
