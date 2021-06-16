NAME         = Chemistry
TARGET       = $${NAME}

QT           = core
QT          += gui
QT          += network
QT          += sql
QT          += script
QT          += positioning
QT          += widgets
QT          += opengl
QT          += printsupport
QT          += multimedia
QT          += multimediawidgets
QT          += opengl
QT          += QtCUDA
QT          += QtOpenCV
QT          += Essentials
QT          += QtCalendar
QT          += QtGMP
QT          += QtGSL
QT          += QtAlgebra
QT          += QtDiscrete
QT          += QtFFT
QT          += Mathematics
QT          += QtFuzzy
QT          += QtFLP
QT          += QtFoundation
QT          += QtGeometry
QT          += QtGadgets
QT          += QtComponents
QT          += QtManagers

load(qt_build_config)
load(qt_module)

INCLUDEPATH += $${PWD}/../../include/$${NAME}

HEADERS     += $${PWD}/../../include/$${NAME}/chemistry.h

SOURCES     += $${PWD}/nPdbViewer.cpp
SOURCES     += $${PWD}/nAtomicPoint.cpp
SOURCES     += $${PWD}/nPdbAtom.cpp
SOURCES     += $${PWD}/nPdbMolecule.cpp

OTHER_FILES += $${PWD}/../../include/$${NAME}/headers.pri

include ($${PWD}/../../doc/Qt/Qt.pri)

TRNAME       = $${NAME}
include ($${PWD}/../../Translations.pri)

win32 {

CONFIG(debug,debug|release) {
  LIBS      += -lfreetyped
} else {
  LIBS      += -lfreetype
}

LIBS        += -lglu32
LIBS        += -lopengl32
}
