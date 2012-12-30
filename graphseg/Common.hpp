/*
 * Common.hpp
 *
 *  Created on: Dec 30, 2012
 *      Author: david
 */

#ifndef DASP_SPECTRAL_COMMON_HPP
#define DASP_SPECTRAL_COMMON_HPP

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
}

#endif
