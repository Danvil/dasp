add_library(pds SHARED
	BlueNoise.cpp
	Grid.cpp
	Mipmaps.cpp
	MultiLayerFloydSteinberg.cpp
	SPDS.cpp
)

include_directories(
)

target_link_libraries(${PROJECT_NAME}
)
