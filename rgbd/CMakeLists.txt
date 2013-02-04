
#
# Project file for rgbd library which provides RGDB streams
#
# parameters:
#  - DASP_HAS_OPENNI: set if OpenNI is available
#

# general settings

SET(rgbd_DEFINES
	-DSLIMAGE_IO_OPENCV
)

SET(rgbd_SOURCES
	rgbd.cpp
)

SET(rgbd_INCLUDES
)

SET(rgbd_LIBRARIES
	opencv_core
	opencv_highgui
	snappy
	boost_filesystem
)


# OpenNI specific settings

if (DASP_HAS_OPENNI)

	LIST(APPEND rgbd_DEFINES
		-DDASP_HAS_OPENNI
	)
	
	LIST(APPEND rgbd_SOURCES
		KinectGrabber.cpp
	)
	
	LIST(APPEND rgbd_INCLUDES
		${OPENNI_INCLUDE_DIR}
	)
	
	LIST(APPEND rgbd_LIBRARIES
		OpenNI
		boost_signals
	)

endif (DASP_HAS_OPENNI)


ADD_DEFINITIONS(${rgbd_DEFINES})
INCLUDE_DIRECTORIES(${rgbd_INCLUDES})
ADD_LIBRARY(rgbd SHARED ${rgbd_SOURCES})
TARGET_LINK_LIBRARIES(rgbd ${rgbd_LIBRARIES})
