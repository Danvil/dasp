add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
)

add_executable(normals_cmd main.cpp)

target_link_libraries(normals_cmd
	libdasp
	opencv_core
	opencv_highgui
	boost_program_options
)

