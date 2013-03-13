add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
)

add_executable(dasp_eval main.cpp)

target_link_libraries(dasp_eval
	libdasp
	opencv_core
	opencv_highgui
	boost_program_options
)

