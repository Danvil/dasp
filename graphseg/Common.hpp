/*
 * Common.hpp
 *
 *  Created on: Dec 30, 2012
 *      Author: david
 */

#ifndef DASP_SPECTRAL_COMMON_HPP
#define DASP_SPECTRAL_COMMON_HPP

#include <boost/graph/adjacency_list.hpp>
#include <Eigen/Dense>

namespace graphseg
{
	namespace detail
	{
		typedef float Real;
		typedef Eigen::MatrixXf Mat;
		typedef Eigen::VectorXf Vec;

		struct EigenComponent
		{
			Real eigenvalue;
			Vec eigenvector;
		};
	}
	
	/** A simple weighted undirected graph */
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
			boost::no_property,
			boost::property<boost::edge_weight_t, float>
	> SpectralGraph;

}

#endif
