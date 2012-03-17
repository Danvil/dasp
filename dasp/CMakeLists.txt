cmake_minimum_required(VERSION 2.6)
link_directories(/home/david/Danvil/lib "$ENV{DANVIL_PATH}"/lib ../libdasp)

project(dasp)

find_package(OpenCV REQUIRED)

file(GLOB my_sources main.cpp)
add_executable(rundasp ${my_sources})
set_target_properties(rundasp PROPERTIES OUTPUT_NAME dasp)

set(CMAKE_CXX_FLAGS "-std=c++0x")

include_directories(
	src
	../../DanvilTools/CT/src
	../../slimage/src
	../libdasp/src
)

target_link_libraries(rundasp
	dasp
	${OpenCV_LIBS}
	boost_thread
	boost_program_options
)

