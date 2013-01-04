add_definitions(-DSLIMAGE_IO_OPENCV)

add_library(dasv SHARED
	dasv.cpp
)

target_link_libraries(dasv
	libdasp
	boost_thread
	opencv_core
	opencv_highgui
)

include_directories(
	${DT_DIR}/CT/src
	${dasp_SOURCE_DIR}/libdasp
	${dasp_SOURCE_DIR}/rgbd
)
