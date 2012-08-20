cmake_minimum_required(VERSION 2.6)

option(DASP_HAS_SIMPLEENGINE "Use DanvilSimpleEngine to enable 3D rendering" OFF) 
option(DASP_HAS_OPENNI "Use OpenNI for Kinect live mode" OFF)

set(OPENNI_INCLUDE_DIR "/home/david/Programs/RGBD/OpenNI/Include" CACHE STRING "OpenNI include directory")
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3" CACHE STRING "Eigen3 include directory")
set(SLIMAGE_INCLUDE_DIR "slimage" CACHE STRING "Slimage include directory")
set(DANVILTOOLS_DIR "/home/david/git/DanvilTools" CACHE STRING "Danvil CT main directory")

add_definitions(-std=c++0x -DBOOST_DISABLE_ASSERTS)

if (DASP_HAS_SIMPLEENGINE)
	link_directories(${DANVILTOOLS_DIR}/SimpleEngine/LinuxRelease)
endif (DASP_HAS_SIMPLEENGINE)

include_directories(
	${EIGEN3_INCLUDE_DIR}
	.
	${DANVILTOOLS_DIR}/CT/src
	${SLIMAGE_INCLUDE_DIR}
)

project(dasp)

add_subdirectory(libdasp)
add_subdirectory(dasp_gui)
add_subdirectory(dasp_tracker)

if (DASP_HAS_OPENNI)
	add_subdirectory(kinect)
endif (DASP_HAS_OPENNI)
