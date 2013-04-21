add_library(pds SHARED
	DeltaDensitySampling.cpp
	Density.cpp
	Fattal.cpp
	FloydSteinberg.cpp
	FloydSteinbergMultiLayer.cpp
	Grid.cpp
	Mipmaps.cpp
	SimplifiedPDS.cpp
	SimplifiedPDSOld.cpp
)

include_directories(
)

target_link_libraries(${PROJECT_NAME}
)
