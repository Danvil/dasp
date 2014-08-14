include_directories(
)

add_executable(density_cmd main.cpp)

target_link_libraries(density_cmd
	density
	opencv_core
	opencv_highgui
	boost_program_options
	boost_timer
	boost_system
)

