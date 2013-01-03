add_library(graphseg SHARED
	Spectral.cpp
	Labeling.cpp
	spectral/SolveSparseTemplate.cpp
)

target_link_libraries(graphseg
	parpack
	arpack
	superlu
)
