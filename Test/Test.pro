#-------------------------------------------------
#
# Project created by QtCreator 2018-04-25T15:51:59
#
#-------------------------------------------------

QT       += widgets testlib

TARGET = tst_testtest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += \
    testmodule.cpp
# install
target.path = $$PWD/tests
INSTALLS += target
