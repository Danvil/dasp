/*
 * Spectral.hpp
 *
 *  Created on: May 19, 2012
 *      Author: david
 */

#ifndef DASP_IMPL_SPECTRAL_HPP_
#define DASP_IMPL_SPECTRAL_HPP_

#include "Spectral/Types.hpp"

namespace dasp
{
	/** Applies spectral graph theory fu to a weighted undirected graph */
	SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev, bool use_dense_solver=false);
}

#endif
