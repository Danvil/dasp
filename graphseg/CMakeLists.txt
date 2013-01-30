
#add_definitions(-DSPECTRAL_VERBOSE)
#add_definitions(-DSEGS_DBG_PRINT)

add_library(graphseg SHARED
	graphseg.cpp
	Labeling.cpp
	spectral/arpackpp.cpp
	spectral/eigen.cpp
	spectral/ietl.cpp
	spectral/lapack.cpp
	spectral/magma.cpp
)

include_directories(graphseg
	/home/david/Downloads/ietl-2.2
)

if (USE_MAGMA)
	add_definitions(-DADD_)
	include_directories(graphseg
		/home/david/Downloads/magma-1.3.0/include
		/usr/local/cuda-5.0/include
	)
endif (USE_MAGMA)

if (SOLVER_IETL)
	include_directories(graphseg
		/home/david/Downloads/ietl-2.2
	)
endif (SOLVER_IETL)

target_link_libraries(graphseg
	parpack
	arpack
	superlu
	lapack
	blas
)

if (USE_MAGMA)
	target_link_libraries(graphseg
		magma
		magmablas
		lapack
		f2c
		cublas
		m
		cudart
		blas
	)
endif (USE_MAGMA)
