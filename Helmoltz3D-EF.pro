#-------------------------------------------------
#
# Project created by QtCreator 2018-02-07T13:33:55
#
#-------------------------------------------------

QT       += core gui 3dcore 3drender 3dinput 3dextras datavisualization

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets


TARGET = Helmoltz3D-EF
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    maillage3d.cpp \
    utils.cpp \
    map_mainwindow.cpp \
    rhs.cpp \
#SOURCES FOR MATH PARSER
    muparserx/parser/mpVariable.cpp \
    muparserx/parser/mpValueCache.cpp \
    muparserx/parser/mpValue.cpp \
    muparserx/parser/mpValReader.cpp \
    muparserx/parser/mpTokenReader.cpp \
    muparserx/parser/mpTest.cpp \
    muparserx/parser/mpScriptTokens.cpp \
    muparserx/parser/mpRPN.cpp \
    muparserx/parser/mpParserMessageProvider.cpp \
    muparserx/parser/mpParserBase.cpp \
    muparserx/parser/mpParser.cpp \
    muparserx/parser/mpPackageUnit.cpp \
    muparserx/parser/mpPackageStr.cpp \
    muparserx/parser/mpPackageNonCmplx.cpp \
    muparserx/parser/mpPackageMatrix.cpp \
    muparserx/parser/mpPackageCommon.cpp \
    muparserx/parser/mpPackageCmplx.cpp \
    muparserx/parser/mpOprtPostfixCommon.cpp \
    muparserx/parser/mpOprtNonCmplx.cpp \
    muparserx/parser/mpOprtMatrix.cpp \
    muparserx/parser/mpOprtIndex.cpp \
    muparserx/parser/mpOprtCmplx.cpp \
    muparserx/parser/mpOprtBinCommon.cpp \
    muparserx/parser/mpOprtBinAssign.cpp \
    muparserx/parser/mpIValue.cpp \
    muparserx/parser/mpIValReader.cpp \
    muparserx/parser/mpIToken.cpp \
    muparserx/parser/mpIPackage.cpp \
    muparserx/parser/mpIOprt.cpp \
    muparserx/parser/mpIfThenElse.cpp \
    muparserx/parser/mpICallback.cpp \
    muparserx/parser/mpFuncStr.cpp \
    muparserx/parser/mpFuncNonCmplx.cpp \
    muparserx/parser/mpFuncMatrix.cpp \
    muparserx/parser/mpFuncCommon.cpp \
    muparserx/parser/mpFuncCmplx.cpp \
    muparserx/parser/mpError.cpp \
    bc.cpp \
    feproblem.cpp \
    MatSparseC3.cpp \
    solver.cpp \
    orbittransformcontroller.cpp \
    viewer3d.cpp \
    scenemodifier.cpp \
    trianglemeshrenderer.cpp \
    scatterdatamodifier.cpp \
    parametersdialog.cpp

HEADERS += \
    maillage3d.h \
    utils.h \
    map_mainwindow.h \
    C3.h \
    rhs.h \
#SOURCES FOR MATH PARSER
    muparserx/parser/utGeneric.h \
    muparserx/parser/suStringTokens.h \
    muparserx/parser/suSortPred.h \
    muparserx/parser/mpVariable.h \
    muparserx/parser/mpValueCache.h \
    muparserx/parser/mpValue.h \
    muparserx/parser/mpValReader.h \
    muparserx/parser/mpTypes.h \
    muparserx/parser/mpTokenReader.h \
    muparserx/parser/mpTest.h \
    muparserx/parser/mpStack.h \
    muparserx/parser/mpScriptTokens.h \
    muparserx/parser/mpRPN.h \
    muparserx/parser/mpParserMessageProvider.h \
    muparserx/parser/mpParserBase.h \
    muparserx/parser/mpParser.h \
    muparserx/parser/mpPackageUnit.h \
    muparserx/parser/mpPackageStr.h \
    muparserx/parser/mpPackageNonCmplx.h \
    muparserx/parser/mpPackageMatrix.h \
    muparserx/parser/mpPackageCommon.h \
    muparserx/parser/mpPackageCmplx.h \
    muparserx/parser/mpOprtPostfixCommon.h \
    muparserx/parser/mpOprtNonCmplx.h \
    muparserx/parser/mpOprtMatrix.h \
    muparserx/parser/mpOprtIndex.h \
    muparserx/parser/mpOprtCmplx.h \
    muparserx/parser/mpOprtBinCommon.h \
    muparserx/parser/mpOprtBinAssign.h \
    muparserx/parser/mpMatrixError.h \
    muparserx/parser/mpMatrix.h \
    muparserx/parser/mpIValue.h \
    muparserx/parser/mpIValReader.h \
    muparserx/parser/mpIToken.h \
    muparserx/parser/mpIPrecedence.h \
    muparserx/parser/mpIPackage.h \
    muparserx/parser/mpIOprt.h \
    muparserx/parser/mpIfThenElse.h \
    muparserx/parser/mpICallback.h \
    muparserx/parser/mpFwdDecl.h \
    muparserx/parser/mpFuncStr.h \
    muparserx/parser/mpFuncNonCmplx.h \
    muparserx/parser/mpFuncMatrix.h \
    muparserx/parser/mpFuncCommon.h \
    muparserx/parser/mpFuncCmplx.h \
    muparserx/parser/mpError.h \
    muparserx/parser/mpDefines.h \
    muparserx/parser/mpCompat.h \
    bc.h \
    feproblem.h \
    RNM/RNM.hpp \
    RNM/RNM_tpl.hpp \
    RNM/RNM_opc.hpp \
    RNM/RNM_op.hpp \
    R3.h \
    R2.h \
    MatSparseC3.h \
    solver.h \
    orbittransformcontroller.h \
    viewer3d.h \
    scenemodifier.h \
    trianglemeshrenderer.h \
    scatterdatamodifier.h \
    parametersdialog.h

CONFIG += c++11

FORMS += \
    map_mainwindow.ui \
    parametersdialog.ui

