add_library(graphseg SHARED
	Spectral.cpp
	Labeling.cpp
	spectral/SolveSparseTemplate.cpp
)

include_directories(graphseg
	/home/david/Downloads/ietl-2.2
)

#link_directories(/opt/acml5.3.0/gfortran64/lib)

target_link_libraries(graphseg
	parpack
	arpack
	superlu
	lapack
	blas
#	acml
)

