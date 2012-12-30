/*
 * Spectral.cpp
 *
 *  Created on: May 19, 2012
 *      Author: david
 */

#include "Spectral.hpp"
#include "spectral/SpectralGraphAnalysis.hpp"

namespace graphseg
{
	SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev, bool use_dense_solver)
	{
		return SpectralGraphAnalysis(graph, num_ev, use_dense_solver);
	}
}
