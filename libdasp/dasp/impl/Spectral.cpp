/*
 * Spectral.cpp
 *
 *  Created on: May 19, 2012
 *      Author: david
 */

#include "Spectral.hpp"
#include "Spectral/SolveTemplate.hpp"
#include <Eigen/Eigenvalues>
#include <boost/assert.hpp>
#include <iostream>

namespace dasp
{

	SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev, bool use_dense_solver)
	{
		return detail::SolveTemplate(graph, num_ev, use_dense_solver);
	}

}
