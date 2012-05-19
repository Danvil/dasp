/*
 * Spectral.hpp
 *
 *  Created on: May 19, 2012
 *      Author: david
 */

#ifndef DASP_IMPL_SPECTRAL_HPP_
#define DASP_IMPL_SPECTRAL_HPP_

#include "../Graph.hpp"
#include <Eigen/Dense>

namespace dasp
{
	/** A simple weighted undirected graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			boost::no_property,
			boost::property<boost::edge_weight_t, float>
	> SpectralGraph;

	/** Applies spectral graph theory fu to a weighted undirected graph */
	SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev);

}

#endif
