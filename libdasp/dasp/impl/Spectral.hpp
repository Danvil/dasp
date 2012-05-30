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
#include <vector>

namespace dasp
{
	/** A simple weighted undirected graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			boost::no_property,
			boost::property<boost::edge_weight_t, float>
	> SpectralGraph;

	/** Applies spectral graph theory fu to a weighted undirected graph */
	SpectralGraph SolveSpectral(const SpectralGraph& graph, unsigned int num_ev, bool use_dense_solver=false);

	namespace detail
	{

		typedef double Real;
		typedef Eigen::MatrixXf Mat;
		typedef Eigen::VectorXf Vec;

		struct EigenComponent
		{
			Real eigenvalue;
			Vec eigenvector;
		};

		typedef std::vector<EigenComponent> PartialEigenSolution;

		PartialEigenSolution SpectralEigenSolveDense(const SpectralGraph& graph, unsigned int num_ev);

		PartialEigenSolution SpectralEigenSolveSparse(const SpectralGraph& graph, unsigned int num_ev);

	}

}

#endif
