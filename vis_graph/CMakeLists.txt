add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
	${DT_DIR}/CT/src
	${CANDY_DIR}/libcandy/src
	${dasp_SOURCE_DIR}/libdasp
	${dasp_SOURCE_DIR}/dasv
	${dasp_SOURCE_DIR}/rgbd
)

add_executable(vis_graph main.cpp)

target_link_libraries(vis_graph
	dasv
	candy
	opencv_core
	opencv_highgui
	boost_program_options
)
