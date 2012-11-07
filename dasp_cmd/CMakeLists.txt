add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
)

add_executable(dasp_cmd main.cpp)

target_link_libraries(dasp_cmd
	libdasp
	opencv_core
	opencv_highgui
	boost_program_options
)

