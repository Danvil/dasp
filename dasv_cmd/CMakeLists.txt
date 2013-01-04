add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
	${dasp_SOURCE_DIR}/dasv
	${dasp_SOURCE_DIR}/rgbd
)

add_executable(dasv_cmd main.cpp)

target_link_libraries(dasv_cmd
	dasv
	rgbd
	opencv_core
	opencv_highgui
	boost_program_options
)

