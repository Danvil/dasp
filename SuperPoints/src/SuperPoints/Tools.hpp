/*
 * Tools.hpp
 *
 *  Created on: Feb 19, 2012
 *      Author: david
 */

#ifndef POINTSANDNORMALS_HPP_
#define POINTSANDNORMALS_HPP_

#include <Slimage/Slimage.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>
#include <ctype.h>

namespace dasp {

/***
 *
 */
struct Camera
{
	float cx, cy;
	float focal;
	float z_slope;

	/** Projects a 3D point into the image plane */
	Eigen::Vector2f project(const Eigen::Vector3f& p) const {
		return Eigen::Vector2f(p[0] / p[2] * focal + cx, p[1] / p[2] * focal + cy);
	}

	/** Computes a 3D point from pixel position and depth */
	Eigen::Vector3f unproject(unsigned int x, unsigned int y, uint16_t depth) const {
		float z = float(depth) * z_slope;
		return Eigen::Vector3f(z*(float(x) - cx)/focal, z*(float(y) - cy)/focal, z);
	}

	/** Gets kinect depth for a 3D point */
	uint16_t depth(const Eigen::Vector3f& p) const {
		return static_cast<uint16_t>(p[2] / z_slope);
	}

};

template<typename K>
inline K LocalFiniteDifferences(K v0, K v1, K v2, K v3, K v4)
{
	K a = std::abs(v2 + v0 - static_cast<K>(2)*v1);
	K b = std::abs(v4 + v2 - static_cast<K>(2)*v3);
	if(a < b) {
		return v2 - v0;
	}
	else {
		return v4 - v2;
	}
}

inline Eigen::Vector2f LocalDepthGradient(const slimage::Image1ui16& depth, unsigned int j, unsigned int i, float base_scale, const Camera& camera)
{
	uint16_t d00 = depth(j,i);

	float h = camera.focal / float(d00);

	// compute w = base_scale*f/d
	unsigned int w = std::max(static_cast<unsigned int>(std::round(base_scale*h)), 2u);

	// can not compute the gradient at the border, so return 0
	if(i < w || depth.height() - w <= i || j < w || depth.width() - w <= j) {
		return Eigen::Vector2f::Zero();
	}

	int dx = LocalFiniteDifferences<int>(
		depth(j-w,i),
		depth(j-w/2,i),
		d00,
		depth(j+w/2,i),
		depth(j+w,i)
	);

	int dy = LocalFiniteDifferences<int>(
		depth(j,i-w),
		depth(j,i-w/2),
		d00,
		depth(j,i+w/2),
		depth(j,i+w)
	);

	// Theoretically scale == base_scale, but w must be an integer, so we
	// compute scale from the actually used w.

	// compute 1 / (scale) = 1 / (w*d/f) = f / (w*d)
	float scl = h / float(w);

	return Eigen::Vector2f(float(dx) * scl, float(dy) * scl);
}

/** Fits a plane into points and returns the plane normal */
template<typename T, typename F>
Eigen::Vector3f FitNormal(const std::vector<T>& points, F f)
{
//		return Eigen::Vector3f(0.0f, 0.0f, 1.0f);
	// compute covariance matrix
//		Eigen::Matrix3f A = Eigen::Matrix3f::Zero();
//		for(const Eigen::Vector3f& p : points) {
//			A += p * p.transpose();
//		}
	float xx=0.0f, xy=0.0f, xz=0.0f, yy=0.0f, yz=0.0f, zz=0.0f;
	for(auto it=points.begin(); it!=points.end(); ++it) {
		const Eigen::Vector3f& p = f(*it);
		float x = p[0];
		float y = p[1];
		float z = p[2];
		xx += x*x;
		xy += x*y;
		xz += x*z;
		yy += y*y;
		yz += y*z;
		zz += z*z;
	};
	Eigen::Matrix3f A; A << xx, xy, xz, xy, yy, yz, xz, yz, zz;
	// compute eigenvalues/-vectors
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver;
	solver.computeDirect(A);
	// take eigenvector (first eigenvalue is smallest!)
	return solver.eigenvectors().col(0).normalized();
}

}

#endif
