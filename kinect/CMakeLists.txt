add_library(kinect SHARED KinectGrabber.cpp)

include_directories(
	${OPENNI_INCLUDE_DIR}
)

target_link_libraries(kinect
	OpenNI
	boost_signals
	boost_thread
)
