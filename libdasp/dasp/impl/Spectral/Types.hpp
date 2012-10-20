/*
 * Types.hpp
 *
 *  Created on: Okt 20, 2012
 *      Author: david
 */

#ifndef DASP_IMPL_SPECTRAL_TYPES_HPP_
#define DASP_IMPL_SPECTRAL_TYPES_HPP_

#include <Eigen/Dense>
#include <boost/graph/adjacency_list.hpp>
#include <vector>

namespace dasp
{
	/** A simple weighted undirected graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			boost::no_property,
			boost::property<boost::edge_weight_t, float>
	> SpectralGraph;

	namespace detail
	{
		extern bool cVerbose;

		typedef double Real;
		typedef Eigen::MatrixXf Mat;
		typedef Eigen::VectorXf Vec;

		struct EigenComponent
		{
			Real eigenvalue;
			Vec eigenvector;
		};

		typedef std::vector<EigenComponent> PartialEigenSolution;
	}

}

#endif
