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
		float z = convertKinectToMeter(depth);
		return Eigen::Vector3f(z*(float(x) - cx)/focal, z*(float(y) - cy)/focal, z);
	}

	/** Gets kinect depth for a 3D point */
	uint16_t depth(const Eigen::Vector3f& p) const {
		return convertMeterToKinect(p[2]);
	}

	/** Convert kinect depth to meter */
	float convertKinectToMeter(uint16_t d) const {
		return static_cast<float>(d) * z_slope;
	}

	float convertKinectToMeter(int d) const {
		return static_cast<float>(d) * z_slope;
	}

	float convertKinectToMeter(float d) const {
		return d * z_slope;
	}

	/** Convert meter to kinect depth */
	uint16_t convertMeterToKinect(float z) const {
		return static_cast<uint16_t>(z / z_slope);
	}

};

template<typename K>
inline float LocalFiniteDifferencesKinect(K v0, K v1, K v2, K v3, K v4)
{
	if(v0 == 0 && v4 == 0 && v1 != 0 && v3 != 0) {
		return float(v3 - v1);
	}

	bool left_invalid = (v0 == 0 || v1 == 0);
	bool right_invalid = (v3 == 0 || v4 == 0);
	if(left_invalid && right_invalid) {
		return 0.0f;
	}
	else if(left_invalid) {
		return float(v4 - v2);
	}
	else if(right_invalid) {
		return float(v2 - v0);
	}
	else {
		float a = static_cast<float>(std::abs(v2 + v0 - static_cast<K>(2)*v1));
		float b = static_cast<float>(std::abs(v4 + v2 - static_cast<K>(2)*v3));
		float p, q;
		if(a + b == 0.0f) {
			p = q = 0.5f;
		}
		else {
			p = a / (a + b);
			q = b / (a + b);
		}
		return q * static_cast<float>(v2 - v0) + p * static_cast<float>(v4 - v2);
	}
}

inline Eigen::Vector2f LocalDepthGradient(const slimage::Image1ui16& depth, unsigned int j, unsigned int i, float base_scale_m, const Camera& camera)
{
	uint16_t d00 = depth(j,i);

	if(d00 == 0) {
		return Eigen::Vector2f::Zero();
	}

	float h = camera.focal / camera.convertKinectToMeter(d00);

	// compute w = base_scale*f/d
	unsigned int w = std::max(static_cast<unsigned int>(std::round(base_scale_m*h)), 4u);
	if(w % 2 == 1) w++;

	// can not compute the gradient at the border, so return 0
	if(i < w || depth.height() - w <= i || j < w || depth.width() - w <= j) {
		return Eigen::Vector2f::Zero();
	}

	float dx = LocalFiniteDifferencesKinect<int>(
		depth(j-w,i),
		depth(j-w/2,i),
		d00,
		depth(j+w/2,i),
		depth(j+w,i)
	);

	float dy = LocalFiniteDifferencesKinect<int>(
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

	Eigen::Vector2f g(camera.convertKinectToMeter(dx) * scl, camera.convertKinectToMeter(dy) * scl);

	const float cGMax = 3.0f;
	g[0] = std::min(+cGMax, std::max(-cGMax, g[0]));
	g[1] = std::min(+cGMax, std::max(-cGMax, g[1]));

	return g;
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
