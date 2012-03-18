cmake_minimum_required(VERSION 2.6)

link_directories(/home/david/Danvil/lib "$ENV{DANVIL_PATH}"/lib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(Eigen3 REQUIRED)

add_definitions(-std=c++0x)

include_directories(
	${EIGEN3_INCLUDE_DIR}
	${OpenCV_INCLUDE_DIRS}
	../DanvilTools/CT/src
	../slimage/src
)

project(dasp)
add_subdirectory(libdasp)
add_subdirectory(dasp)
add_subdirectory(boundary_recall)
add_subdirectory(redcarpet)

