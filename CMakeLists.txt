cmake_minimum_required(VERSION 2.6)

option(DASP_HAS_CANDY "Use DanvilSimpleEngine to enable 3D rendering" OFF) 
option(DASP_HAS_OPENNI "Use OpenNI for Kinect live mode" OFF)

set(OPENNI_INCLUDE_DIR "/home/david/Programs/RGBD/OpenNI/Include" CACHE STRING "OpenNI include directory")
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3" CACHE STRING "Eigen3 include directory")
set(SLIMAGE_INCLUDE_DIR "slimage" CACHE STRING "Slimage include directory")
set(DT_DIR "/home/david/git/DanvilTools" CACHE STRING "Danvil tools main directory")
set(CANDY_DIR "/home/david/git/candy" CACHE STRING "Danvil candy 3D engine main directory")

add_definitions(-std=c++0x -DBOOST_DISABLE_ASSERTS)

if (DASP_HAS_CANDY)
	link_directories(/home/david/Projects/candy/release/libcandy) # FIXME
endif (DASP_HAS_CANDY)

include_directories(
	${EIGEN3_INCLUDE_DIR}
	.
	${SLIMAGE_INCLUDE_DIR}
)

project(dasp)

add_subdirectory(libdasp)
add_subdirectory(dasp_cmd)
add_subdirectory(dasp_gui)

add_subdirectory(dasv)
add_subdirectory(dasv_cmd)

if (DASP_HAS_OPENNI)
	add_subdirectory(kinect)
endif (DASP_HAS_OPENNI)
