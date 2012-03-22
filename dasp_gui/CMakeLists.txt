find_package(Qt4 REQUIRED)

add_definitions(-DSLIMAGE_IO_QT)

set(OPENNI_INCLUDE_DIR /home/david/Programs/RGBD/OpenNI/Include)
set(OPENNI_LIB_DIR D:/Program Files/OpenNI/Lib)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
	${OPENNI_INCLUDE_DIR}
)

SET(QT_USE_QTOPENGL TRUE)

SET(dasp_gui_SOURCES
	main.cpp
	WdgtSuperpixelParameters.cpp
	WdgtKinectSuperPoints.cpp
)
if (DASP_HAS_SIMPLEENGINE)
SET(dasp_gui_HEADERS
	WdgtKinectSuperPoints.h
	WdgtSuperpixelParameters.h
	../../DanvilTools/SimpleEngine/src/Danvil/SimpleEngine/System/GLSystemQtWindow.h
)
else (DASP_HAS_SIMPLEENGINE)
SET(dasp_gui_HEADERS
	WdgtKinectSuperPoints.h
	WdgtSuperpixelParameters.h
)
endif (DASP_HAS_SIMPLEENGINE)

SET(dasp_gui_FORMS
	WdgtSuperpixelParameters.ui
	WdgtKinectSuperPoints.ui
)

QT4_WRAP_CPP(dasp_gui_HEADERS_MOC ${dasp_gui_HEADERS})
QT4_WRAP_UI(dasp_gui_FORMS_HEADERS ${dasp_gui_FORMS})

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS} -DQT_NO_KEYWORDS)

ADD_EXECUTABLE(dasp_gui
	${dasp_gui_SOURCES} ${dasp_gui_HEADERS_MOC} ${dasp_gui_FORMS_HEADERS}
	KinectGrabber.cpp
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TARGET_LINK_LIBRARIES(dasp_gui
	libdasp
	${QT_LIBRARIES}
	OpenNI
	boost_signals
	boost_thread
)

if (DASP_HAS_SIMPLEENGINE)
	add_definitions(-DDASP_HAS_SIMPLEENGINE)
	include_directories(../../DanvilTools/SimpleEngine/src)
	target_link_libraries(dasp_gui DanvilSimpleEngine)
endif (DASP_HAS_SIMPLEENGINE)

