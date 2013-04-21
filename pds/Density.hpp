#ifndef INCLUDED_PDS_DENSITY_HPP
#define INCLUDED_PDS_DENSITY_HPP

#include <Eigen/Dense>
#include <vector>

namespace pds
{
	Eigen::MatrixXf PointDensity(const std::vector<Eigen::Vector2f>& seeds, const Eigen::MatrixXf& target);
}

#endif
