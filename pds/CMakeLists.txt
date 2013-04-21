add_library(pds SHARED
	BlueNoise.cpp
	FloydSteinberg.cpp
	Grid.cpp
	Mipmaps.cpp
	MultiLayerFloydSteinberg.cpp
	SimplifiedPoissonDiskSampling.cpp
	SimplifiedPoissonDiskSamplingOld.cpp
)

include_directories(
)

target_link_libraries(${PROJECT_NAME}
)
