/*
 * PointsAndNormals.hpp
 *
 *  Created on: Feb 2, 2012
 *      Author: david
 */

#ifndef POINTSANDNORMALS_HPP_
#define POINTSANDNORMALS_HPP_

#include <Slimage/Slimage.hpp>
#include <Slimage/Parallel.h>
#include <eigen3/Eigen/Dense>
#include <list>
#include <vector>

namespace dasp {
namespace PointsAndNormals {

constexpr float cZSlope = 0.001f;
constexpr float cFocal = 580.0f;
constexpr float cWindowMeters = 0.04f;

/** Computes a 3D point from pixel position and depth */
inline
Eigen::Vector3f PointFromDepth(unsigned int x, unsigned int y, uint16_t depth, unsigned int width, unsigned int height) {
	float z = float(depth) * cZSlope;
	return Eigen::Vector3f(
			z*(float(x) - float(width/2))/cFocal,
			z*(float(y) - float(height/2))/cFocal,
			z);
}

inline
Eigen::Vector3f PointFromZ(unsigned int x, unsigned int y, float z, unsigned int width, unsigned int height) {
	return Eigen::Vector3f(
			z*(float(x) - float(width/2))/cFocal,
			z*(float(y) - float(height/2))/cFocal,
			z);
}

/** Computes 3D points for every pixel from depth values */
slimage::Image3f ComputePoints(const slimage::Image1ui16& depth, slimage::ThreadingOptions opt);

namespace NormalModes {
	enum NormalMode {
		KNearestEigen,
		EightMeanCross
	};
}
typedef NormalModes::NormalMode NormalMode;

/** Computes normals for every pixel from 3D points (and depth) */
slimage::Image3f ComputeNormals(const slimage::Image1ui16& depth, const slimage::Image3f& points, slimage::ThreadingOptions opt, NormalMode mode=NormalModes::KNearestEigen);

}}

#endif
