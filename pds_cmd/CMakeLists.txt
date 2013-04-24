ADD_DEFINITIONS(-DSLIMAGE_IO_OPENCV)

include_directories(
)

add_executable(pds_cmd main.cpp)

target_link_libraries(pds_cmd
	pds
	opencv_core
	opencv_highgui
	boost_program_options
	boost_timer
	boost_system
)

