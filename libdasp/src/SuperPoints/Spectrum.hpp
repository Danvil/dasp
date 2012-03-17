/*
 * Spectrum.hpp
 *
 *  Created on: Mar 2, 2012
 *      Author: david
 */

#ifndef SPECTRUM_HPP_
#define SPECTRUM_HPP_
//----------------------------------------------------------------------------//
#include <eigen3/Eigen/Dense>
#include <vector>
//----------------------------------------------------------------------------//
namespace Dasp {
//----------------------------------------------------------------------------//

namespace Spectrum
{

	std::vector<float> PointsToDistances(const std::vector<Eigen::Vector3f>& points);

};

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
