#ifndef PLOTTING_HPP_
#define PLOTTING_HPP_
//----------------------------------------------------------------------------//
#include <Danvil/Color.h>
#include <Eigen/Dense>
//----------------------------------------------------------------------------//
namespace common {
//----------------------------------------------------------------------------//

/** Computes a color to express a similarity value between 0 and 1 */
inline Eigen::Vector3f SimilarityColor(float x)
{
	static auto cm = Danvil::ContinuousIntervalColorMapping<float, float>::Factor_Black_Blue_Red_Yellow_White();
	cm.setRange(0.0f, 1.0f);
	Danvil::Colorf color = cm(x);
	return {color.r,color.g,color.b};
}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
