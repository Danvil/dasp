add_definitions(-DSLIMAGE_IO_OPENCV)

include_directories(
)

add_library(rgbd SHARED rgbd.cpp)

target_link_libraries(rgbd
	kinect
	opencv_core
	opencv_highgui
)

