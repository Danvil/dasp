include_directories(
	${dasp_SOURCE_DIR}/libdasp
)

add_executable(dasp_density main.cpp)

target_link_libraries(dasp_density
	libdasp
	opencv_core
	opencv_highgui
	boost_program_options
	boost_system
)
