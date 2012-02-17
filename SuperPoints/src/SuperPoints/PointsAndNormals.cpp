/*
 * PointsAndNormals.cpp
 *
 *  Created on: Feb 8, 2012
 *      Author: david
 */

#include "PointsAndNormals.hpp"
#include <Slimage/Slimage.hpp>
#include <Slimage/Parallel.h>
#include <eigen3/Eigen/Dense>
#include <list>
#include <vector>
//----------------------------------------------------------------------------//
namespace dasp {
namespace PointsAndNormals {
//----------------------------------------------------------------------------//

inline
Eigen::Vector3f ReadVector3fFromMem(const float* p_points) {
	return Eigen::Vector3f(p_points[0], p_points[1], p_points[2]);
}

inline
void WriteVector3fToMem(float* p_points, const Eigen::Vector3f& pnt) {
	p_points[0] = pnt[0];
	p_points[1] = pnt[1];
	p_points[2] = pnt[2];
}

slimage::Image3f ComputePoints(const slimage::Image1ui16& depth, slimage::ThreadingOptions opt)
{
	const unsigned int width = depth.width();
	const unsigned int height = depth.height();
	slimage::Image3f points(width, height);

	const uint16_t* p_depth_begin = depth.begin();
	slimage::ParallelProcess(depth, points, [p_depth_begin,width,height](const uint16_t* p_depth, float* p_points) {
		unsigned int i = p_depth - p_depth_begin;
		unsigned int x = i % width;
		unsigned int y = i / width;
		Eigen::Vector3f pnt = PointFromDepth(x, y, *p_depth, width, height);
		WriteVector3fToMem(p_points, pnt);
	}, opt);

	return points;
}

/** Fits a plane into points and returns the plane normal */
template<typename C>
Eigen::Vector3f FitNormal(const C& points) {
//		return Eigen::Vector3f(0.0f, 0.0f, 1.0f);
	// compute covariance matrix
//		Eigen::Matrix3f A = Eigen::Matrix3f::Zero();
//		for(const Eigen::Vector3f& p : points) {
//			A += p * p.transpose();
//		}
	float xx=0.0f, xy=0.0f, xz=0.0f, yy=0.0f, yz=0.0f, zz=0.0f;
	for(const Eigen::Vector3f& p : points) {
		float x = p[0];
		float y = p[1];
		float z = p[2];
		xx += x*x;
		xy += x*y;
		xz += x*z;
		yy += y*y;
		yz += y*z;
		zz += z*z;
	}
	Eigen::Matrix3f A; A << xx, xy, xz, xy, yy, yz, xz, yz, zz;
	// compute eigenvalues/-vectors
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver;
	solver.computeDirect(A);
	// take eigenvector (first eigenvalue is smallest!)
	return solver.eigenvectors().col(0).normalized();
}

inline
bool InRange(int x, int y, unsigned int width, unsigned int height) {
	return (0 <= x && x < int(width) && 0 <= y && y < int(height));
}

inline
bool InRange(unsigned int x, unsigned int y, unsigned int width, unsigned int height) {
	return (x < width && y < height);
}

slimage::Image3f ComputeNormals_KNearestEigen(const slimage::Image1ui16& depth, const slimage::Image3f& points, slimage::ThreadingOptions opt)
{
	struct Vec3fWithNorm {
		Vec3fWithNorm() {}
		Vec3fWithNorm(const Eigen::Vector3f& x)
		: x_(x), norm_2_(x.squaredNorm()) {}
		Eigen::Vector3f x_;
		float norm_2_;
	};

	const unsigned int width = depth.width();
	const unsigned int height = depth.height();
//	const unsigned int K = 80;
//	const int R = 5; // checks (2*R+1)^2 = 121 points
	slimage::Image3f normals(width, height);

	const float* p_point_begin = points.begin();
	slimage::ParallelProcess(points, normals, [&depth,p_point_begin,width,height](const float* p_point, float* p_normal) {
//			WriteVector3fToMem(p_normal, Eigen::Vector3f(0.0f, 0.0f, -1.0f));

		int i = (p_point - p_point_begin)/3;
		int x0 = i % width;
		int y0 = i / width;

		float s = cWindowMeters * cFocal / (float(depth[i]) * cZSlope);
		int R = int(s * 0.5f);
		R = std::max(1, R);
		const float dcut = 0.58742364f * cWindowMeters; // magic number!
		const float dcut2 = dcut * dcut;

		// do not process border
		if(x0-R < 0 || int(width) <= x0+R || y0-R < 0 || int(height) <= y0+R) {
			WriteVector3fToMem(p_normal, Eigen::Vector3f::Zero());
			return;
		}
		// look in neighbourhood in spiral
		Eigen::Vector3f center = ReadVector3fFromMem(p_point);
		if(center.squaredNorm() < 0.01f) {
			WriteVector3fToMem(p_normal, Eigen::Vector3f::Zero());
			return;
		}
		// gather all points in a window
		std::vector<Vec3fWithNorm> selected_points;
		selected_points.reserve((2*R + 1)*(2*R + 1));
		for(int y=-R; y<=+R; y++) {
			for(int x=-R; x<=+R; x++) {
				Eigen::Vector3f p = ReadVector3fFromMem(p_point_begin + 3*(x0 + x + (y0 + y)*width));
				if(p.squaredNorm() > 0.01f) {
					Vec3fWithNorm dpc(p - center);
					if(dpc.norm_2_ < dcut2) {
						selected_points.push_back(dpc);
					}
				}
			}
		}
		// find k nearest points
#if 0
		unsigned int k = std::min<unsigned int>(K, selected_points.size());
		std::nth_element(selected_points.begin(), selected_points.begin() + k, selected_points.end(),
				[](const Vec3fWithNorm& x, const Vec3fWithNorm& y) { return x.norm_2_ < y.norm_2_; }
		);
#else
		unsigned int k = selected_points.size();
#endif
		if(k < 5) {
			WriteVector3fToMem(p_normal, Eigen::Vector3f::Zero());
			return;
		}
		std::vector<Eigen::Vector3f> k_nearest(k);
		auto it = selected_points.begin();
		for(unsigned int i=0; i<k; i++, ++it) {
				k_nearest[i] = it->x_;
		}
		// fit plane into points
		Eigen::Vector3f normal = FitNormal(k_nearest);
		// normal shall face to camera
		if(normal[2] > 0) {
			normal = -normal;
		}
		// write normal
		WriteVector3fToMem(p_normal, normal);
	}, opt);

	return normals;
}

slimage::Image4f ComputeIntegralImage(const slimage::Image3f& points)
{
	unsigned int w = points.width() + 1;
	unsigned int h = points.height() + 1;
	slimage::Image4f integral(w, h);
	for(unsigned int x=0; x<w; x++) {
		integral(x,0) = slimage::Pixel4f{{0.0f,0.0f,0.0f,0.0f}};
	}
	for(unsigned int y=1; y<h; y++) {
		slimage::Pixel4f row_sum{{0.0f,0.0f,0.0f,0.0f}};
		integral(0,y) = row_sum;
		for(unsigned int x=1; x<w; x++) {
			auto v = points(x-1,y-1);
			row_sum[0] += v[0];
			row_sum[1] += v[1];
			row_sum[2] += v[2];
			bool is_zero = (v[0] == 0.0f && v[1] == 0.0f && v[2] == 0.0f);
			row_sum[3] += (is_zero ? 0.0f : 1.0f);
			integral(x,y) = row_sum + integral(x,y-1);
		}
	}
	return integral;
}

Eigen::Vector4f ToEigen(const slimage::detail::PixelAccess<float,4>& x) {
	return Eigen::Vector4f(x[0], x[1], x[2], x[3]);
}

Eigen::Vector3f ComputeNormal(const Eigen::Vector4f& u, const Eigen::Vector4f& v, const Eigen::Vector3f& center, unsigned int& cnt) {
	if(u[3] == 0.0f || v[3] == 0.0f) {
		return Eigen::Vector3f::Zero();
	}
	cnt ++;
	Eigen::Vector3f un(u[0]/u[3], u[1]/u[3], u[2]/u[3]);
	Eigen::Vector3f vn(v[0]/v[3], v[1]/v[3], v[2]/v[3]);
	return (un - center).cross(vn - center).normalized();
}

slimage::Image3f ComputeNormals_EightMeanCross(const slimage::Image1ui16& depth, const slimage::Image3f& points, slimage::ThreadingOptions opt)
{
	slimage::Image4f integral = ComputeIntegralImage(points);

	const unsigned int width = depth.width();
	const unsigned int height = depth.height();
	slimage::Image3f normals(width, height);

	const float* p_normal_begin = normals.begin();
	slimage::ParallelProcess(depth, normals, [p_normal_begin,width,height,&integral](const uint16_t* p_depth, float* p_normal) {
		int i = (p_normal - p_normal_begin)/3;
		int x0 = i % width;
		int y0 = i / width;
		// compute depth adaptive scale
		float s = cWindowMeters * cFocal / (float(*p_depth) * cZSlope);
		int R = int(s * 0.5f);
		R = std::max(1, R);
		const int p[4] = {-3*R-1, -R, +R+1, +3*R+2 };
		if(x0+p[0] < 0 || int(width + 1) <= x0+p[3] || y0+p[0] < 0 || int(height + 1) <= y0+p[3]) {
			WriteVector3fToMem(p_normal, Eigen::Vector3f::Zero());
			return;
		}
		const Eigen::Vector4f I[4][4] = {
				{ ToEigen(integral(x0+p[0], y0+p[0])), ToEigen(integral(x0+p[1], y0+p[0])), ToEigen(integral(x0+p[2], y0+p[0])), ToEigen(integral(x0+p[3], y0+p[0])) },
				{ ToEigen(integral(x0+p[0], y0+p[1])), ToEigen(integral(x0+p[1], y0+p[1])), ToEigen(integral(x0+p[2], y0+p[1])), ToEigen(integral(x0+p[3], y0+p[1])) },
				{ ToEigen(integral(x0+p[0], y0+p[2])), ToEigen(integral(x0+p[1], y0+p[2])), ToEigen(integral(x0+p[2], y0+p[2])), ToEigen(integral(x0+p[3], y0+p[2])) },
				{ ToEigen(integral(x0+p[0], y0+p[3])), ToEigen(integral(x0+p[1], y0+p[3])), ToEigen(integral(x0+p[2], y0+p[3])), ToEigen(integral(x0+p[3], y0+p[3])) }
		};
		Eigen::Vector4f M[3][3];
		for(unsigned int y=0; y<3; y++) {
			for(unsigned int x=0; x<3; x++) {
				M[y][x] = I[y+1][x+1] + I[y][x] - I[y][x+1] - I[y+1][x];
			}
		}
		// compute mean position of center area
		Eigen::Vector4f center4 = M[1][1];
		if(center4[3] == 0.0f) {
			WriteVector3fToMem(p_normal, Eigen::Vector3f::Zero());
			return;
		}
		Eigen::Vector3f center(center4[0]/center4[3], center4[1]/center4[3], center4[2]/center4[3]);
		// compute normal as mean of some cross products
		Eigen::Vector3f n = Eigen::Vector3f::Zero();
		unsigned int cnt = 0;
		n += ComputeNormal(M[0][2], M[0][0], center, cnt);
		n += ComputeNormal(M[2][0], M[2][2], center, cnt);
		n += ComputeNormal(M[0][1], M[1][0], center, cnt);
		n += ComputeNormal(M[2][1], M[1][2], center, cnt);
		if(cnt == 0) {
			WriteVector3fToMem(p_normal, Eigen::Vector3f::Zero());
			return;
		}
		// normalize normal
		WriteVector3fToMem(p_normal, n.normalized());
//		WriteVector3fToMem(p_normal, Eigen::Vector3f(0.0f, 0.0f, -1.0f));
	}, opt);

	return normals;
}

slimage::Image3f ComputeNormals(const slimage::Image1ui16& depth, const slimage::Image3f& points, slimage::ThreadingOptions opt, NormalMode mode)
{
	switch(mode) {
	case NormalModes::KNearestEigen:
		return ComputeNormals_KNearestEigen(depth, points, opt);
	case NormalModes::EightMeanCross:
		return ComputeNormals_EightMeanCross(depth, points, opt);
	default:
		// TODO invalid mode
		return slimage::Image3f();
	}
}

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
