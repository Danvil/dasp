/*
 * Spectral.hpp
 *
 *  Created on: Dec 30, 2012
 *      Author: david
 */

#ifndef DASP_SPCTRAL_HPP
#define DASP_SPCTRAL_HPP

#include <boost/graph/adjacency_list.hpp>

namespace graphseg
{

	namespace detail
	{
		extern bool cVerbose;
	}

	/** A simple weighted undirected graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			boost::no_property,
			boost::property<boost::edge_weight_t, float>
	> SpectralGraph;

	/** Applies spectral graph theory fu to a weighted, undirected graph */
	SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev, bool use_dense_solver=false);

	/** Applies MCL graph segmentation to a weighted, undirected graph */
	SpectralGraph SolveMCL(const SpectralGraph& graph, float p, unsigned int iterations);

}

#endif
