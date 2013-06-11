ADD_DEFINITIONS(-DSLIMAGE_IO_OPENCV)

add_library(density SHARED
	PointDensity.cpp
	ScalePyramid.cpp
	Smooth.cpp
	Visualization.cpp
)

target_link_libraries(density
	opencv_core
	opencv_highgui
)

