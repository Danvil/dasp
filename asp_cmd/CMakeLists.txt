ADD_DEFINITIONS(-DSLIMAGE_IO_OPENCV)

include_directories(
)

add_executable(asp_cmd main.cpp)

target_link_libraries(asp_cmd
	asp
	pds
	opencv_core
	opencv_highgui
	boost_program_options
	boost_timer
	boost_system
)

