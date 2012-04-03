cmake_minimum_required(VERSION 2.6)

link_directories(/home/david/Danvil/lib "$ENV{DANVIL_PATH}"/lib)

option(DASP_HAS_SIMPLEENGINE "Use DanvilSimpleEngine to enable 3D rendering" ON) 

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /home/david/git/SuperPoints/cmake/Modules/)
#find_package(Eigen3 REQUIRED)

add_definitions(-std=c++0x -DBOOST_DISABLE_ASSERTS -O3)

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


include_directories(
#	${EIGEN3_INCLUDE_DIR}
	/usr/local/include/eigen3
	${OpenCV_INCLUDE_DIRS}
	../DanvilTools/CT/src
	../slimage/src
)

project(dasp)
add_subdirectory(libdasp)
add_subdirectory(dasp_gui)
add_subdirectory(dasp_cmd)
add_subdirectory(eval_superpixels)
add_subdirectory(eval_dasp_segments)
add_subdirectory(eval_gPb)
add_subdirectory(boundary_to_labels)
add_subdirectory(redcarpet)
