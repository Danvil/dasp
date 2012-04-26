/*
 * PointSet.hpp
 *
 *  Created on: Apr 20, 2012
 *      Author: david
 */

#ifndef POINTSET_HPP_
#define POINTSET_HPP_

#include <Eigen/Dense>
#include <vector>

struct Point
{
	Eigen::Vector3f position;
	Eigen::Vector3f normal;
	Eigen::Vector3f color;

	friend Point operator*(const Eigen::Affine3f& T, const Point& s) {
		return {
			T * s.position, // transform position
			T.rotation() * s.normal,  // transform normal only with rotation!
			s.color };
	}
};

typedef std::vector<Point> PointSet;

inline float Distance(const Point& p, const Point& q) {
//		const float cSuperpixelRadius = 0.02f;
//		return (p.position - q.position).squaredNorm() / (4.0f * cSuperpixelRadius * cSuperpixelRadius) + (p.color - q.color).squaredNorm();
	return (p.position - q.position).squaredNorm();
}

inline std::size_t FindClosestPoint(const PointSet& points, const Point& p) {
	float best_distance = 1e9;
	std::size_t best_id = 0;
	for(std::size_t i=0; i<points.size(); i++) {
		float d = Distance(p, points[i]);
		if(d < best_distance) {
			best_distance = d;
			best_id = i;
		}
	}
	return best_id;
}

inline PointSet operator*(const Eigen::Affine3f& T, const PointSet& s) {
	PointSet u(s.size());
	for(unsigned int i=0; i<s.size(); i++) {
		u[i] = T * s[i];
	}
	return u;
}

#endif
