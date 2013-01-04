find_package(Qt4 REQUIRED)

add_definitions(-DSLIMAGE_IO_QT)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
	${dasp_SOURCE_DIR}/dasv
	${dasp_SOURCE_DIR}/rgbd
	${DT_DIR}/CT/src
	${CANDY_DIR}/libcandy/src
)

SET(QT_USE_QTOPENGL TRUE)

SET(dasv_vis_SOURCES
	main.cpp
	WdgtDasvVis.cpp
)
SET(dasv_vis_HEADERS
	WdgtDasvVis.h
	${CANDY_DIR}/libcandy/src/Candy/System/GLSystemQtWindow.h
)

SET(dasv_vis_FORMS
	WdgtDasvVis.ui
)

QT4_WRAP_CPP(dasv_vis_HEADERS_MOC ${dasv_vis_HEADERS})
QT4_WRAP_UI(dasv_vis_FORMS_HEADERS ${dasv_vis_FORMS})

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS} -DQT_NO_KEYWORDS)

ADD_EXECUTABLE(dasv_vis
	${dasv_vis_SOURCES} ${dasv_vis_HEADERS_MOC} ${dasv_vis_FORMS_HEADERS}
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TARGET_LINK_LIBRARIES(dasv_vis
	dasv
	rgbd
	candy
	${QT_LIBRARIES}
	boost_program_options
)
