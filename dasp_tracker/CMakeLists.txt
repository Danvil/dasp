add_library(dasp_tracker SHARED
	dasp/DaspTracker.cpp
)

include_directories(
	${dasp_SOURCE_DIR}/libdasp
)

target_link_libraries(dasp_tracker
	libdasp
)

if (DASP_HAS_SIMPLEENGINE)
	add_definitions(-DDASP_HAS_SIMPLEENGINE)
	include_directories(${DANVILTOOLS_DIR}/SimpleEngine/src)
	target_link_libraries(dasp_tracker DanvilSimpleEngine)
endif (DASP_HAS_SIMPLEENGINE)

# add opencv (used in Segmentation.cpp)
add_definitions(-DSLIMAGE_IO_OPENCV)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(dasp_tracker ${OpenCV_LIBS})
