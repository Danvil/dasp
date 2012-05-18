/*
 * Metric.hpp
 *
 *  Created on: May 18, 2012
 *      Author: david
 */

#ifndef DASP_METRIC_HPP_
#define DASP_METRIC_HPP_

#include "Point.hpp"
#include <Eigen/Dense>

namespace dasp
{

	namespace metric
	{
		inline float ImageDistanceRaw(const Point& x, const Point& y) {
			return (x.pos - y.pos).squaredNorm() / (y.image_super_radius * y.image_super_radius);
		}

		inline float SpatialDistanceRaw(const Point& x, const Point& y) {
			return (x.world - y.world).squaredNorm();// / (y.spatial_normalizer * y.spatial_normalizer);
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

//	struct ClassicNeighbourhoodMetric
//	{
//		float operator()(const Point& x, const Point& y) {
//			NeighbourhoodGraphEdgeData& edge = neighbourhood_graph[eid];
//			// compute cost using color and normal
//			edge.c_px = metric::ImageDistanceRaw(x, y);
//			edge.c_world = metric::SpatialDistanceRaw(x, y) / (opt.base_radius * opt.base_radius); // F IXME HAAAACK
//			edge.c_color = metric::ColorDistanceRaw(x, y);
//			edge.c_normal = metric::NormalDistanceRaw(x, y);
//			// F IXME metric needs central place!
//			if(opt.weight_image == 0.0f) {
//				// dasp
//				MetricDASP fnc(opt.weight_spatial, opt.weight_color, opt.weight_normal, opt.base_radius);
//				edge.weight = fnc(x, y);
//			}
//			else {
//				// slic
//				MetricSLIC fnc(opt.weight_image, opt.weight_color);
//				edge.weight = fnc(x, y);
//			}
//		}
//	};

	template<bool SupressConvexEdges=true>
	struct ClassicSpectralAffinity
	{
		ClassicSpectralAffinity(unsigned int num_superpixels, float superpixel_radius, float w_spatial=1.0f, float w_color=1.0f, float w_normal=1.0f)
		: num_superpixels_(num_superpixels),
		  superpixel_radius_(superpixel_radius),
		  w_spatial(w_spatial),
		  w_color(w_color),
		  w_normal(w_normal)
		{}

		float operator()(const Point& x, const Point& y) const {
			float c_world = metric::SpatialDistanceRaw(x, y) / (superpixel_radius_ * superpixel_radius_); // FIXME HAAAACK
			float c_color = metric::ColorDistanceRaw(x, y);
			float w_maha_color = c_color / (std::sqrt(static_cast<float>(num_superpixels_)) * cWeightRho);
			float w_maha_spatial= std::max(0.0f, std::min(4.0f, c_world/4.0f - 1.0f)); // distance of 2 indicates normal distance
			float w_maha_normal;
			if(SupressConvexEdges) {
				// only use concave edges
				Eigen::Vector3f d = (y.world - x.world).normalized();
				float ca1 = x.computeNormal().dot(d);
				float ca2 = y.computeNormal().dot(d);
				float w = ca1 - ca2;
				//float w = std::acos(ca2) - std::acos(ca1);
				if(w < 0.0f) {
					w_maha_normal = 0.0f;
				}
				else {
					w_maha_normal = 3.0f * w;
					//w_maha_normal = 3.0f * (1 - std::cos(w));
				}
			}
			else {
				w_maha_normal = 3.0f*metric::NormalDistanceRaw(x, y);
			}
			// compute total edge connectivity
			return std::exp(-(w_spatial*w_maha_spatial + w_color*w_maha_color + w_normal*w_maha_normal));
		}

	private:
		static constexpr float cWeightRho = 0.01f; // 640x480 clusters would yield 0.1 which is used in gPb
		unsigned int num_superpixels_;
		float superpixel_radius_;
		float w_spatial;
		float w_color;
		float w_normal;
	};

	struct ClassicSpectralAffinitySLIC
	{
		ClassicSpectralAffinitySLIC(unsigned int num_superpixels, float w_color=1.0f)
		: num_superpixels_(num_superpixels),
		  w_color(w_color)
		{}

		float operator()(const Point& x, const Point& y) const {
			float w_maha_color = 4.0f * metric::ColorDistanceRaw(x, y) / (std::sqrt(static_cast<float>(num_superpixels_)) * cWeightRho);
			return std::exp(-w_color*w_maha_color);
		}

	private:
		static constexpr float cWeightRho = 0.01f; // 640x480 clusters would yield 0.1 which is used in gPb
		unsigned int num_superpixels_;
		float w_color;
	};

}

#endif
