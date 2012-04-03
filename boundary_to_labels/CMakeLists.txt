add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
	${OpenCV_INCLUDE_DIRS}
)

find_package(OpenCV REQUIRED)

add_executable(boundary_to_labels boundary_to_labels.cpp)

target_link_libraries(boundary_to_labels
	libdasp
	${OpenCV_LIBS}
)

