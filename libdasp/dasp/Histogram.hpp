/*
 * Histogram.hpp
 *
 *  Created on: Feb 17, 2012
 *      Author: david
 */

#ifndef HISTOGRAM_HPP_
#define HISTOGRAM_HPP_
//----------------------------------------------------------------------------//
#include <eigen3/Eigen/Dense>
#include <cmath>
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

inline
float HistogramDistance(const Eigen::VectorXf& x, const Eigen::VectorXf& y, const Eigen::MatrixXf& A, float m=0.7f)
{
	/// For more details on this distance see the paper:
	///  The Quadratic-Chi Histogram Distance Family
	///  Ofir Pele, Michael Werman
	///  ECCV 2010
	Eigen::VectorXf z = A * (x + y);
	Eigen::VectorXf d = x - y;
	for(int i=0; i<z.rows(); i++) {
		float u = z[i];
		if(u == 0.0f) {
			u = 1.0f;
		}
		else {
			u = std::pow(u, m);
		}
		d[i] /= u;
	}
	float d2 = d.transpose() * A * d;
	return std::sqrt(std::max(0.0f, d2));

}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
