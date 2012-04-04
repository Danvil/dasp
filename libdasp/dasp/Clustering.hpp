/*
 * Clustering.hpp
 *
 *  Created on: Apr 4, 2012
 *      Author: david
 */

#ifndef DASP_CLUSTERING_HPP_
#define DASP_CLUSTERING_HPP_

#include "Point.hpp"
#include "Parameters.hpp"
#include <Slimage/Slimage.hpp>
#include <eigen3/Eigen/Dense>

namespace dasp
{

	namespace metric
	{
		inline float ImageDistanceRaw(const Point& x, const Point& y) {
			return (x.pos - y.pos).squaredNorm() / (y.image_super_radius * y.image_super_radius);
		}

		inline float SpatialDistanceRaw(const Point& x, const Point& y) {
			return (x.world - y.world).squaredNorm();
		}

		inline float ColorDistanceRaw(const Point& u, const Point& v) {
			return (u.color - v.color).squaredNorm();
		}

		inline float NormalDistanceRaw(const Point& u, const Point& v) {
			// we want to compute 1 - dot(n(u), n(v))
			// the normal is implicitly given by the gradient
			// multiplying with the circularity yields the required normalization
			return 1.0f - (u.gradient.dot(v.gradient) + 1.0f) * u.circularity * v.circularity;
		}
	}

	struct MetricDASP
	{
		MetricDASP(float w_r, float w_c, float w_n, float R) {
			weights_[0] = w_r / (R * R);
			weights_[1] = w_c;
			weights_[2] = w_n;
		}

		float operator()(const Point& p, const Point& q) const {
			return weights_.dot(
					Eigen::Vector3f(
							metric::SpatialDistanceRaw(p, q),
							metric::ColorDistanceRaw(p, q),
							metric::NormalDistanceRaw(p, q)));
		}

	private:
		Eigen::Vector3f weights_;
	};

	struct MetricSLIC
	{
		MetricSLIC(float w_s, float w_c) {
			weights_[0] = w_s;
			weights_[1] = w_c;
		}

		float operator()(const Point& p, const Point& q) const {
			return weights_.dot(
					Eigen::Vector2f(
							metric::ImageDistanceRaw(p, q),
							metric::ColorDistanceRaw(p, q)));
		}

	private:
		Eigen::Vector2f weights_;
	};

	template<typename METRIC>
	slimage::Image1f ComputeEdges(const ImagePoints& points, const METRIC& mf)
	{
		const unsigned int width = points.width();
		const unsigned int height = points.height();
		slimage::Image1f edges(width, height, slimage::Pixel1f{1e6f});
		// compute edges strength
		for(unsigned int y=1; y<height-1; y++) {
			for(unsigned int x=1; x<width-1; x++) {
				float v;
				const Point& p0 = points(x, y);
				const Point& px1 = points(x-1, y);
				const Point& px2 = points(x+1, y);
				const Point& py1 = points(x, y-1);
				const Point& py2 = points(x, y+1);
				if(p0.isInvalid() || px1.isInvalid() || px2.isInvalid() || py1.isInvalid() || py2.isInvalid()) {
					v = 1e6; // dont want to be here
				}
				else {
					float dx = mf(px1, px2);
					float dy = mf(py1, py2);
					v = dx + dy;
				}
				edges(x,y) = v;
			}
		}
		return edges;
	}

	template<typename METRIC>
	slimage::Image1i IterateClusters(const std::vector<Cluster>& clusters, const ImagePoints& points, const Parameters& opt, const METRIC& mf)
	{
		slimage::Image1i labels(points.width(), points.height(), slimage::Pixel1i{-1});
		std::vector<float> v_dist(points.size(), 1e9);
		// for each cluster check possible points
		for(unsigned int j=0; j<clusters.size(); j++) {
			const Cluster& c = clusters[j];
			int cx = c.center.spatial_x();
			int cy = c.center.spatial_y();
			const int R = int(c.center.image_super_radius * opt.coverage);
			const unsigned int xmin = std::max(0, cx - R);
			const unsigned int xmax = std::min(int(points.width()-1), cx + R);
			const unsigned int ymin = std::max(0, cy - R);
			const unsigned int ymax = std::min(int(points.height()-1), cy + R);
			//unsigned int pnt_index = points.index(xmin,ymin);
			for(unsigned int y=ymin; y<=ymax; y++) {
				for(unsigned int x=xmin; x<=xmax; x++/*, pnt_index++*/) {
					unsigned int pnt_index = points.index(x, y);
					const Point& p = points[pnt_index];
					if(p.isInvalid()) {
						// omit invalid points
						continue;
					}
					float dist = mf(p, c.center);
					float& v_dist_best = v_dist[pnt_index];
					if(dist < v_dist_best) {
						v_dist_best = dist;
						labels[pnt_index] = j;
					}
				}
				//pnt_index -= (xmax - xmin + 1);
				//pnt_index += points.width();
			}
		}
		return labels;
	}

//	template<bool cUseSqrt>
//	inline float ImageDistance(const Point& x, const Point& y) {
//		float d_img;
//		if(cUseSqrt) {
//			d_img = (x.pos - y.pos).norm();
//			d_img /= y.image_super_radius;
//		}
//		else {
//			d_img = (x.pos - y.pos).squaredNorm();
//			d_img /= y.image_super_radius * y.image_super_radius;
//		}
//		return d_img;
//	}
//
//	template<bool cUseSqrt>
//	inline float SpatialDistance(const Point& x, const Point& y, const Parameters& opt) {
//		float d_world;
//		if(cUseSqrt) {
//			d_world = (x.world - y.world).norm();
//			d_world /= opt.base_radius;
//		}
//		else {
//			d_world = (x.world - y.world).squaredNorm();
//			d_world /= opt.base_radius * opt.base_radius;
//		}
//		return d_world;
//	}
//
//	template<bool cUseSqrt, bool cLimit>
//	inline float ColorDistance(const Point& u, const Point& v) {
//		float d_color;
//		if(cUseSqrt) {
//			d_color = (u.color - v.color).norm();
//			if(cLimit) {
//				d_color /= 0.25f;
//			}
//		}
//		else {
//			d_color = (u.color - v.color).squaredNorm();
//			if(cLimit) {
//				d_color /= (0.25f * 0.25f);
//			}
//		}
//		if(cLimit) {
//			if(d_color > 1.0f) {
//				d_color = 1.0f;
//			}
//		}
//		return d_color;
//	}
//
//	inline float NormalDistance(const Point& u, const Point& v) {
//		return NormalDistanceRaw(u, v);
//	}
//
//	inline float DepthWeight(const Point& u, const Point& v) {
//		return 2.0f / (1.0f + 0.5f * (u.depth() + v.depth()));
//	}
//
//	template<bool cUseSqrt=false, bool WeightNormalsByDepth=false, bool cLimitColorDist=false>
//	inline float Distance(const Point& u, const Point& v, const Parameters& opt) const
//	{
//		float d_image = ImageDistance<cUseSqrt>(u,v);
//		float d_world = SpatialDistance<cUseSqrt>(u,v);
//		float d_color = ColorDistance<cUseSqrt,cLimitColorDist>(u,v);
//		float d_normal = NormalDistance(u, v);
//		if(WeightNormalsByDepth) {
//			d_normal *= DepthWeight(u, v);
//		}
//		return opt.weight_spatial * d_world
//			+ opt.weight_image * d_image
//			+ opt.weight_color * d_color
//			+ opt.weight_normal * d_normal;
//	}
//
//	template<bool WeightNormalsByDepth=false>
//	inline float DistanceSpatialColorNormal(const Cluster& x, const Cluster& y, const Parameters& opt) {
//		// world position of cluster centers.
//		float d_world = SpatialDistance<true>(x.center, y.center);
//		// mean cluster colors,
//		float d_color = ColorDistance<true,true>(x.center, y.center);
//		// mean cluster normals and
//		float d_normal = NormalDistance(x.center, y.center);
//		if(WeightNormalsByDepth) {
//			d_normal *= DepthWeight(x.center, y.center);
//		}
//		return opt.weight_spatial*d_world
//			+ opt.weight_color*d_color
//			+ opt.weight_normal*d_normal;
//	}
//
//	template<bool WeightNormalsByDepth=false>
//	inline float DistanceColorNormal(const Cluster& x, const Cluster& y, const Parameters& opt) {
//		// mean cluster colors,
//		float d_color = ColorDistance<true,true>(x.center, y.center);
//		// mean cluster normals and
//		float d_normal = NormalDistance(x.center, y.center);
//		if(WeightNormalsByDepth) {
//			d_normal *= DepthWeight(x.center, y.center);
//		}
//		return opt.weight_color*d_color
//			+ opt.weight_normal*d_normal;
//	}
//
//	template<bool WeightNormalsByDepth=false>
//	inline float DistanceNormal(const Cluster& x, const Cluster& y, const Parameters& opt) {
//		// mean cluster normals and
//		float d_normal = NormalDistance(x.center, y.center);
//		if(WeightNormalsByDepth) {
//			d_normal *= 2.0f / (1.0f + 0.5f * (x.center.depth() + y.center.depth()));
//		}
//		return opt.weight_normal*d_normal;
//	}

}

#endif
