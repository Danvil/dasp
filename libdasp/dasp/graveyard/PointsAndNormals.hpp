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
#include <eigen3/Eigen/Eigenvalues>
#include <list>
#include <vector>
#include <ctype.h>

namespace dasp {

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

}

#endif
