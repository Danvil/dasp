add_library(libdasp SHARED
	dasp/tools/BlueNoise.cpp
	dasp/tools/Mipmaps.cpp
	dasp/tools/Recall.cpp
	dasp/tools/RepairDepth.cpp
	dasp/DaspTracker.cpp
	dasp/Neighbourhood.cpp
	dasp/Plots.cpp
	dasp/Sampling.cpp
	dasp/Segmentation.cpp
	dasp/Superpixels.cpp
)

set_target_properties(libdasp PROPERTIES OUTPUT_NAME dasp)

target_link_libraries(libdasp
	boost_thread
)

if (DASP_HAS_SIMPLEENGINE)
	add_definitions(-DDASP_HAS_SIMPLEENGINE)
	include_directories(${DANVILTOOLS_DIR}/SimpleEngine/src)
	target_link_libraries(libdasp DanvilSimpleEngine)
endif (DASP_HAS_SIMPLEENGINE)

# add opencv (used in Segmentation.cpp)
add_definitions(-DSLIMAGE_IO_OPENCV)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(libdasp ${OpenCV_LIBS})
