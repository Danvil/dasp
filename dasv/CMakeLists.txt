add_definitions(-DSLIMAGE_IO_OPENCV)

add_library(dasv SHARED
	dasv.cpp
)

target_link_libraries(dasv
	boost_thread
	opencv_core
	opencv_highgui
)

#if (DASP_HAS_CANDY)
#	add_definitions(-DDASP_HAS_CANDY)
#	include_directories(
#		${DT_DIR}/CT/src
#		${CANDY_DIR}/libcandy/src
#	)
#	target_link_libraries(libdasp candy)
#endif (DASP_HAS_CANDY)
