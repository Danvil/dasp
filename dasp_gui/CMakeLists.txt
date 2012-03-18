find_package(Qt4 REQUIRED)

add_definitions(-DSLIMAGE_IO_QT)

include_directories(
	${dasp_SOURCE_DIR}/libdasp/src
)

SET(QT_USE_QTOPENGL TRUE)

SET(dasp_gui_SOURCES main.cpp WdgtSuperpixelParameters.cpp WdgtKinectSuperPoints.cpp)
SET(dasp_gui_HEADERS WdgtKinectSuperPoints.h WdgtSuperpixelParameters.h)
SET(dasp_gui_FORMS WdgtSuperpixelParameters.ui WdgtKinectSuperPoints.ui)

QT4_WRAP_CPP(dasp_gui_HEADERS_MOC ${dasp_gui_HEADERS})
QT4_WRAP_UI(dasp_gui_FORMS_HEADERS ${dasp_gui_FORMS})

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS})

ADD_EXECUTABLE(dasp_gui ${dasp_gui_SOURCES} ${dasp_gui_HEADERS_MOC} ${dasp_gui_FORMS_HEADERS})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TARGET_LINK_LIBRARIES(dasp_gui ${QT_LIBRARIES})

