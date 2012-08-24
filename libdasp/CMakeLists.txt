add_library(libdasp SHARED
	dasp/impl/BlueNoise.cpp
	dasp/impl/Mipmaps.cpp
	dasp/impl/RepairDepth.cpp
	dasp/impl/Spectral.cpp
	dasp/impl/SpectralEigenSolveDense.cpp
	dasp/impl/SpectralEigenSolveSparse.cpp
	dasp/Neighbourhood.cpp
	dasp/Plots.cpp
	dasp/Sampling.cpp
	dasp/Segmentation.cpp
	dasp/Superpixels.cpp
)

set_target_properties(libdasp PROPERTIES OUTPUT_NAME dasp)

target_link_libraries(libdasp
	boost_thread
	parpack
	arpack
	superlu
)

if (DASP_HAS_CANDY)
	add_definitions(-DDASP_HAS_CANDY)
	include_directories(
		${DT_DIR}/CT/src
		${CANDY_DIR}/libcandy/src
	)
	target_link_libraries(libdasp candy)
endif (DASP_HAS_CANDY)

# add opencv (used in Segmentation.cpp)
add_definitions(-DSLIMAGE_IO_OPENCV)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(libdasp ${OpenCV_LIBS})
