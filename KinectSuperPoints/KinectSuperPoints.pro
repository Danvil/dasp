TEMPLATE = app
QT += core \
    gui \
    opengl
CONFIG += no_keywords
GIT_BASE = /home/david/git
HEADERS += src/SuperPoints/KinectSuperPoints/WdgtKinectSuperPoints.h \
    src/SuperPoints/KinectSuperPoints/WdgtSuperpixelParameters.h \
    $$GIT_BASE/DanvilTools/SimpleEngine/src/Danvil/SimpleEngine/System/GLSystemQtWindow.h
SOURCES += src/SuperPoints/KinectSuperPoints/WdgtKinectSuperPoints.cpp \
    src/SuperPoints/KinectSuperPoints/WdgtSuperpixelParameters.cpp \
    src/SuperPoints/KinectSuperPoints/main.cpp
FORMS += src/SuperPoints/KinectSuperPoints/WdgtKinectSuperPoints.ui \
    src/SuperPoints/KinectSuperPoints/WdgtSuperpixelParameters.ui
RESOURCES += 
QMAKE_CXXFLAGS += -std=c++0x
INCLUDEPATH += src \
    $$GIT_BASE/SuperPoints/SuperPoints/src \
    $$GIT_BASE/Romeo/RomeoKinect/src \
    $$GIT_BASE/DanvilTools/Voxello/src \
    $$GIT_BASE/DanvilTools/SimpleEngine/src \
    $$GIT_BASE/DanvilTools/CT/src \
    $$GIT_BASE/slimage/src \
    $$GIT_BASE/DanvilTools/Images/src \
    /home/david/Programs/RGBD/OpenNI/Include \
    /usr/local/include/eigen3
CONFIG(debug, debug|release) { 
    DEFINES += DEBUG
    TARGET = KinectSuperPoints-Debug
    LIBS += -lRomeoKinect-Debug \
        -lDanvilImages-Debug \
        -lDanvilSimpleEngine-Debug \
        -lVoxello-Debug \
        -lSuperPoints-Debug
}
CONFIG(release, debug|release) { 
    DEFINES += NDEBUG
    TARGET = KinectSuperPoints
    LIBS += -lRomeoKinect \
        -lDanvilImages \
        -lDanvilSimpleEngine \
        -lVoxello \
        -lSuperPoints
}
LIBS += -L$$(DANVIL_PATH)/lib
LIBS += -lboost_signals \
    -lboost_program_options \
    -lboost_thread \
    -lOpenNI
