file(GLOB my_sources dasp/*.cpp)
add_library(libdasp SHARED ${my_sources})

set_target_properties(libdasp PROPERTIES OUTPUT_NAME dasp)

target_link_libraries(libdasp
	boost_thread
)

if (DASP_HAS_SIMPLEENGINE)
	add_definitions(-DDASP_HAS_SIMPLEENGINE)
	include_directories(../../DanvilTools/SimpleEngine/src)
	target_link_libraries(libdasp DanvilSimpleEngine)
endif (DASP_HAS_SIMPLEENGINE)

