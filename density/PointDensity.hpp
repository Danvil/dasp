#ifndef INCLUDED_PDS_DENSITY_HPP
#define INCLUDED_PDS_DENSITY_HPP

#include <Eigen/Dense>
#include <vector>

namespace density
{
	Eigen::MatrixXf PointDensity(const std::vector<Eigen::Vector2f>& points, const Eigen::MatrixXf& density);
}

#endif
