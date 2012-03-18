include_directories(
	src
	../../DanvilTools/SimpleEngine/src
)

file(GLOB my_sources src/**/*.cpp)
add_library(libdasp SHARED ${my_sources})

set_target_properties(libdasp PROPERTIES OUTPUT_NAME dasp)

target_link_libraries(libdasp
	DanvilSimpleEngine
)

