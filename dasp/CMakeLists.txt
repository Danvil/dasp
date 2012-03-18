add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
	src
	${dasp_SOURCE_DIR}/libdasp/src
	${OpenCV_INCLUDE_DIRS}
)

find_package(OpenCV REQUIRED)

file(GLOB my_sources main.cpp)
add_executable(rundasp ${my_sources})

set_target_properties(rundasp PROPERTIES OUTPUT_NAME dasp)

target_link_libraries(rundasp
	libdasp
	${OpenCV_LIBS}
	boost_thread
	boost_program_options
)

