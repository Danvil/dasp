

add_library(graphseg SHARED
	Spectral.cpp
	Labeling.cpp
	spectral/SolveSparseTemplate.cpp
	spectral/Magma.cpp
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
