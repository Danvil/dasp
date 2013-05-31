add_library(pds SHARED
	DeltaDensitySampling.cpp
	Fattal.cpp
	FloydSteinberg.cpp
	FloydSteinbergMultiLayer.cpp
	Grid.cpp
	SimplifiedPDS.cpp
	SimplifiedPDSOld.cpp
)

target_link_libraries(pds
	density
)

