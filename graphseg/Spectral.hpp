/*
 * Spectral.hpp
 *
 *  Created on: Dec 30, 2012
 *      Author: david
 */

#ifndef DASP_SPCTRAL_HPP
#define DASP_SPCTRAL_HPP

#include "Common.hpp"

namespace graphseg
{

	namespace detail
	{
		extern bool cVerbose;
	}

	enum class SpectralMethod {
		Eigen,
		Magma,
		ArpackPPSparse,
		Lapack,
		Ietl
	};

	/** Applies spectral graph theory fu to a weighted, undirected graph */
	SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev, SpectralMethod method);

	/** Applies MCL graph segmentation to a weighted, undirected graph */
	SpectralGraph SolveMCL(const SpectralGraph& graph, float p, unsigned int iterations);

}

#endif
