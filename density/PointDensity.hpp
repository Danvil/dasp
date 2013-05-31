#ifndef INCLUDED_PDS_DENSITY_HPP
#define INCLUDED_PDS_DENSITY_HPP

#include <Eigen/Dense>
#include <vector>

namespace density
{

	/** Gaussian blur using density-adaptive radius as kernel radius */
	Eigen::MatrixXf DensityAdaptiveSmooth(const Eigen::MatrixXf& d);

	/** Computes density approximation for a set of points  */
	Eigen::MatrixXf PointDensity(const std::vector<Eigen::Vector2f>& points, const Eigen::MatrixXf& density);

}

#endif
