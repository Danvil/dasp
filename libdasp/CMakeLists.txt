cmake_minimum_required(VERSION 2.6)
link_directories(/home/david/Danvil/lib "$ENV{DANVIL_PATH}"/lib)

project(dasp)

file(GLOB my_sources src/**/*.cpp)
add_library(libdasp SHARED ${my_sources})
set_target_properties(libdasp PROPERTIES OUTPUT_NAME dasp)

set(CMAKE_CXX_FLAGS "-std=c++0x")

include_directories(
	src
	../../DanvilTools/CT/src
	../../slimage/src
	../../DanvilTools/SimpleEngine/src
	/usr/local/include/eigen3
)

target_link_libraries(libdasp
	DanvilSimpleEngine
)

