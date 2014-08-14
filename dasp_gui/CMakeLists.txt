find_package(Qt4 REQUIRED)

add_definitions(-DSLIMAGE_IO_QT)

if (DASP_HAS_CANDY)
	add_definitions(-DDASP_HAS_CANDY)
endif (DASP_HAS_CANDY)

if (DASP_HAS_OPENNI)
	add_definitions(-DDASP_HAS_OPENNI)
endif (DASP_HAS_OPENNI)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
	${dasp_SOURCE_DIR}/rgbd
)

SET(QT_USE_QTOPENGL TRUE)

SET(dasp_gui_SOURCES
	main.cpp
	DaspProcessing.cpp
	WdgtKinectSuperPoints.cpp
	WdgtSettings.cpp
	../common/WdgtDaspParameters.cpp
	WdgtBenchmark.cpp
	WdgtAbout.cpp
)
if (DASP_HAS_CANDY)
	SET(dasp_gui_HEADERS
		WdgtKinectSuperPoints.h
		WdgtSettings.h
		../common/WdgtDaspParameters.h
		WdgtBenchmark.h
		WdgtAbout.h
		${CANDY_DIR}/libcandy/Candy/System/GLSystemQtWindow.h
	)
else (DASP_HAS_CANDY)
	SET(dasp_gui_HEADERS
		WdgtKinectSuperPoints.h
		WdgtSettings.h
		../common/WdgtDaspParameters.h
		WdgtBenchmark.h
		WdgtAbout.h
	)
endif (DASP_HAS_CANDY)

SET(dasp_gui_FORMS
	WdgtKinectSuperPoints.ui
	WdgtSettings.ui
	../common/WdgtDaspParameters.ui
	WdgtBenchmark.ui
	WdgtAbout.ui
)

QT4_WRAP_CPP(dasp_gui_HEADERS_MOC ${dasp_gui_HEADERS})
QT4_WRAP_UI(dasp_gui_FORMS_HEADERS ${dasp_gui_FORMS})

INCLUDE(${QT_USE_FILE})
ADD_DEFINITIONS(${QT_DEFINITIONS} -DQT_NO_KEYWORDS)

ADD_EXECUTABLE(dasp_gui
	${dasp_gui_SOURCES} ${dasp_gui_HEADERS_MOC} ${dasp_gui_FORMS_HEADERS}
)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

TARGET_LINK_LIBRARIES(dasp_gui
	libdasp
	rgbd
	${QT_LIBRARIES}
	boost_signals
	boost_thread
	boost_program_options
	boost_system
	pthread
)

if (DASP_HAS_CANDY)
	add_definitions(-DDASP_HAS_CANDY)
	include_directories(
		${CANDY_DIR}/libcandy
	)
	target_link_libraries(dasp_gui candy)
endif (DASP_HAS_CANDY)
