add_definitions(-DSLIMAGE_IO_OPENCV)
add_definitions(-DDASP_HAS_SIMPLEENGINE)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
	${dasp_SOURCE_DIR}/kinect
	${OPENNI_INCLUDE_DIR}
	../../DanvilTools/SimpleEngine/src
	${OpenCV_INCLUDE_DIRS}
)

find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

add_executable(register register.cpp DaspRegistration.cpp
	../../DanvilTools/SimpleEngine/src/Danvil/SimpleEngine/System/GlutSystem.cpp
)

target_link_libraries(register
	libdasp
	kinect
	DanvilSimpleEngine
	${OpenCV_LIBS}
	boost_program_options
	${GLUT_LIBRARY}
	${OPENGL_LIBRARY}
	m
)
