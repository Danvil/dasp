/*
 * PointsAndNormals.hpp
 *
 *  Created on: Feb 2, 2012
 *      Author: david
 */

#ifndef POINTSANDNORMALS_HPP_
#define POINTSANDNORMALS_HPP_

#include <Danvil/Images/Image.h>
#include <Danvil/Images/ImageOps.h>
#include <Danvil/Images/Parallel.h>
#include <eigen3/Eigen/Dense>
#include <list>
#include <vector>

namespace PointsAndNormals
{
	inline
	Eigen::Vector3f PointFromDepth(unsigned int x, unsigned int y, uint16_t depth, unsigned int width, unsigned int height) {
		const float z_slope = 0.001f;
		const float scl = 1.0f / 580.0f;
		float z = float(depth) * z_slope;
		return Eigen::Vector3f(
				z*(float(x) - float(width/2))*scl,
				z*(float(y) - float(height/2))*scl,
				z);
	}

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
		return (0 <= x && x < width && 0 <= y && y < height);
	}

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

	void ComputePoints(const Danvil::Images::Image1ui16Ptr& depth, Danvil::Images::Image3fPtr& points, Danvil::Images::ThreadingOptions opt) {
		const unsigned int width = points->width();
		const unsigned int height = points->height();
		Danvil::Images::ImageOps::Resize(points, depth);
		const uint16_t* p_depth_begin = depth->begin();
		Danvil::Images::ParallelProcess(depth, points, [p_depth_begin,width,height](const uint16_t* p_depth, float* p_points) {
			unsigned int i = p_depth - p_depth_begin;
			unsigned int x = i % width;
			unsigned int y = i / width;
			Eigen::Vector3f pnt = PointFromDepth(x, y, *p_depth, width, height);
			WriteVector3fToMem(p_points, pnt);
		}, opt);
	}

	struct Vec3fWithNorm {
		Vec3fWithNorm() {}
		Vec3fWithNorm(const Eigen::Vector3f& x)
		: x_(x), norm_2_(x.squaredNorm()) {}
		Eigen::Vector3f x_;
		float norm_2_;
	};

	inline
	void InsertSorted(std::list<Vec3fWithNorm>& v, const Eigen::Vector3f& p) {
		Vec3fWithNorm pwn(p);
		auto it = std::lower_bound(v.begin(), v.end(), pwn, [](const Vec3fWithNorm& x, const Vec3fWithNorm& y) {
			return x.norm_2_ < y.norm_2_;
		});
		v.insert(it, pwn);
	}

	inline
	void Insert(std::list<Vec3fWithNorm>& v, const float* p_point, const Eigen::Vector3f& center, int x, int y, unsigned int width, unsigned int height) {
		Eigen::Vector3f p = ReadVector3fFromMem(p_point + 3*(x + y*width));
		if(p.squaredNorm() > 0.01f) {
			InsertSorted(v, p - center);
		}
	}

	inline
	void InsertIfValid(std::list<Vec3fWithNorm>& v, const float* p_point, const Eigen::Vector3f& center, int x, int y, unsigned int width, unsigned int height) {
		if(InRange(x, y, width, height)) {
			Insert(v, p_point, center, x, y, width, height);
		}
	}

	void ComputeNormals(const Danvil::Images::Image1ui16Ptr& depth, const Danvil::Images::Image3fPtr& points, Danvil::Images::Image3fPtr& normals, Danvil::Images::ThreadingOptions opt) {
		const unsigned int width = points->width();
		const unsigned int height = points->height();
		const unsigned int K = 32;
		const int R = 3; // checks (2*R+1)^2 = 121 points
		Danvil::Images::ImageOps::Resize(normals, points);
		const float* p_point_begin = points->begin();
		Danvil::Images::ParallelProcess(points, normals, [p_point_begin,width,height](const float* p_point, float* p_normal) {
//			WriteVector3fToMem(p_normal, Eigen::Vector3f(0.0f, 0.0f, -1.0f));

			int i = (p_point - p_point_begin)/3;
			int x0 = i % width;
			int y0 = i / width;
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
			std::list<Vec3fWithNorm> sorted_points;
			sorted_points.push_back(Vec3fWithNorm(Eigen::Vector3f::Zero()));
			for(int r=1; r<=R; r++) {
				for(int x=-r; x<=+r; x++) {
					// insert (x0+x, y0-r) and (x0+x, y0+r)
					Insert(sorted_points, p_point_begin, center, x0+x, y0-r, width, height);
					Insert(sorted_points, p_point_begin, center, x0+x, y0+r, width, height);
				}
				for(int y=-r+1; y<=+r-1; y++) {
					// insert (x0-r, y0+y) and (x0+r, y0+y)
					Insert(sorted_points, p_point_begin, center, x0-r, y0+y, width, height);
					Insert(sorted_points, p_point_begin, center, x0+r, y0+y, width, height);
				}
			}
			if(sorted_points.size() < 3) {
				WriteVector3fToMem(p_normal, Eigen::Vector3f::Zero());
				return;
			}
			// fit plane into points
			unsigned int k = std::min<unsigned int>(K, sorted_points.size());
			std::vector<Eigen::Vector3f> k_nearest(k);
			auto it = sorted_points.begin();
			for(unsigned int i=0; i<k; i++, ++it) {
					k_nearest[i] = it->x_;
			}
			Eigen::Vector3f normal = FitNormal(k_nearest);
			// normal shall face to camera
			if(normal[2] > 0) {
				normal = -normal;
			}
			WriteVector3fToMem(p_normal, normal);
		}, opt);
	}
}

#endif /* POINTSANDNORMALS_HPP_ */
