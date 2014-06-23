
SET(graphseg_DEFINES
#	-DSPECTRAL_VERBOSE
#	-DSEGS_DBG_PRINT
)

if (USE_SOLVER_ARPACK)
	LIST(APPEND graphseg_DEFINES
		-DUSE_SOLVER_ARPACK
	)
endif (USE_SOLVER_ARPACK)

if (USE_SOLVER_MAGMA)
	LIST(APPEND graphseg_DEFINES
		-DUSE_SOLVER_MAGMA
		-DADD_
	)
endif (USE_SOLVER_MAGMA)

if (USE_SOLVER_IETL)
	LIST(APPEND graphseg_DEFINES
		-DUSE_USE_SOLVER_IETL
	)
endif (USE_SOLVER_IETL)

ADD_DEFINITIONS(${graphseg_DEFINES})


SET(graphseg_SOURCES
	graphseg.cpp
	Labeling.cpp
	spectral/eigen.cpp
	spectral/lapack.cpp
)

if (USE_SOLVER_ARPACK)
	LIST(APPEND graphseg_SOURCES spectral/arpackpp.cpp)
endif (USE_SOLVER_ARPACK)

if (USE_SOLVER_MAGMA)
	LIST(APPEND graphseg_SOURCES spectral/magma.cpp)
endif (USE_SOLVER_MAGMA)

if (USE_SOLVER_IETL)
	LIST(APPEND graphseg_SOURCES spectral/ietl.cpp)
endif (USE_SOLVER_IETL)

ADD_LIBRARY(graphseg SHARED ${graphseg_SOURCES})


SET(graphseg_INCLUDES
)

if (USE_SOLVER_MAGMA)
	LIST(APPEND graphseg_INCLUDES
		/home/david/Downloads/magma-1.3.0/include
		/usr/local/cuda-5.0/include
	)
endif (USE_SOLVER_MAGMA)

if (USE_SOLVER_IETL)
	LIST(APPEND graphseg_INCLUDES
		/home/david/Downloads/ietl-2.2
	)
endif (USE_SOLVER_IETL)

INCLUDE_DIRECTORIES(graphseg ${graphseg_INCLUDES})


SET(graphseg_LIBRARIES
	lapack
	blas
)

if (USE_SOLVER_ARPACK)
	LIST(APPEND graphseg_LIBRARIES
		parpack
		arpack
		superlu
	)
endif (USE_SOLVER_ARPACK)

if (USE_SOLVER_MAGMA)
	LIST(APPEND graphseg_LIBRARIES
		magma
		magmablas
		lapack
		f2c
		cublas
		m
		cudart
		blas
	)
endif (USE_SOLVER_MAGMA)

TARGET_LINK_LIBRARIES(graphseg ${graphseg_LIBRARIES})
