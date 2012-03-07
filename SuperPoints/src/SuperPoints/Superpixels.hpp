/*
 * Superpixels.hpp
 *
 *  Created on: Jan 26, 2012
 *      Author: david
 */

#ifndef SUPERPIXELS_HPP_
#define SUPERPIXELS_HPP_

#include "Tools.hpp"
#include "SuperpixelHistogram.hpp"
#include <Slimage/Slimage.hpp>
#include <Slimage/Parallel.h>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>

namespace dasp
{
	struct Seed
	{
		int x, y;
		float scala;
	};

	struct Point
	{
		/** pixel image position */
		Eigen::Vector2f pos;

		/** point color */
		Eigen::Vector3f color;

		/** depth in mm as integer (0 if invalid) */
		uint16_t depth_i16;

		/** position [m] of world source point */
		Eigen::Vector3f world;

		/** local depth gradient (depth[m]/distance[m]) */
		Eigen::Vector2f gradient;

		/** direction of local normal */
		Eigen::Vector3f normal;

		/** estimated radius [px] on the image screen of a super pixel at point depth */
		float image_super_radius;

		/** Invalid points have a kinect depth of 0 */
		bool isInvalid() const {
			return depth_i16 == 0;
		}

		/** Valid points have a kinect depth > 0 */
		bool isValid() const {
			return depth_i16 > 0;
		}

		/** Image x coordinate [px] */
		int spatial_x() const {
			return static_cast<int>(pos[0] + 0.5f);
		}

		/** Image y coordinate [px] */
		int spatial_y() const {
			return static_cast<int>(pos[1] + 0.5f);
		}

		/** Depth [m] of point */
		float depth() const {
			return world[2];
		}

	public:
		 EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	};

	struct ImagePoints
	{
		ImagePoints()
		: width_(0), height_(0) {}
		ImagePoints(unsigned int width, unsigned int height)
		: width_(width), height_(height), points_(width*height) {
			auto it = points_.begin();
			for(unsigned int y=0; y<height; y++) {
				for(unsigned int x=0; x<width; x++, ++it) {
					Point& p = *it;
					p.pos[0] = float(x);
					p.pos[1] = float(y);
				}
			}
		}
		unsigned int width() const {
			return width_;
		}
		unsigned int height() const {
			return height_;
		}
		unsigned int size() const {
			return points_.size();
		}
		size_t index(unsigned int x, unsigned int y) const {
			return x + y*width_;
		}
		size_t index(const Point& p) const {
			return index(p.spatial_x(), p.spatial_y());
		}
		const Point& operator()(const Eigen::Vector2f& p) const {
			int x = static_cast<int>(p[0] + 0.5f);
			int y = static_cast<int>(p[1] + 0.5f);
			return (*this)(x,y);
		}
		const Point& operator()(unsigned int x, unsigned int y) const {
			return points_[index(x,y)];
		}
		Point& operator()(unsigned int x, unsigned int y) {
			return points_[index(x,y)];
		}
		const Point& operator[](unsigned int i) const {
			return points_[i];
		}
		Point& operator[](unsigned int i) {
			return points_[i];
		}
	private:
		unsigned int width_, height_;
		std::vector<Point> points_;
	};

	namespace SeedModes {
		enum Type {
			EquiDistant,
			DepthShooting,
			DepthMipmap,
			DepthBlueNoise,
			DepthFloyd,
			Delta
		};
	}
	typedef SeedModes::Type SeedMode;

	struct Parameters
	{
		Parameters() {
			weight_color = 1.0f;
			weight_spatial = 1.0f;
			weight_normal = 1.0f;
			weight_depth = 1.0f;
			iterations = 3;
			coverage = 1.7f;
			base_radius = 0.02f;
			seed_mode = SeedModes::DepthMipmap;
			gradient_adaptive_density = true;
		}

		/** camera parameters */
		Camera camera;

		float weight_color;
		float weight_spatial;
		float weight_normal;
		float weight_depth;

		/** Number of iterations for superpixel k-means clustering */
		unsigned int iterations;

		/** Superpixel cluster search radius factor */
		float coverage;

		/** Desired radius of a surface element */
		float base_radius;

		/** Method used to compute seed points */
		SeedMode seed_mode;

		bool gradient_adaptive_density;

		/** Pixel scala at depth
		 * Radius [px] of a surface element of size base radius [m] and
		 * at given depth [kinect] on the image sensor
		 */
		float computePixelScala(uint16_t depth) const {
			return (depth == 0) ? 0.0f : (camera.focal / camera.convertKinectToMeter(depth) * base_radius);
		}
	};

	struct Cluster
	{
		Point center;

		std::vector<unsigned int> pixel_ids;

		bool hasPoints() const {
			return pixel_ids.size() > 3;
		}

		/** Eigenvalues of the covariance matrix */
		float t, b, a;
		/** eccentricity of the ellipse described by a and b */
		float eccentricity;
		/** pi*a*b */
		float area;
		/** sqrt(a*b) */
		float radius;
		/** number of pixel which are within superpixel radius but are not part of the superpixel */
		float coverage;

		void UpdateCenter(const ImagePoints& points, const Camera& cam);

		void ComputeClusterInfo(const ImagePoints& points, const Parameters& opt);

	};

	struct ClusterGroupInfo
	{
		Histogram<float> hist_eccentricity;
		Histogram<float> hist_radius;
		Histogram<float> hist_thickness;
	};

	class Clustering
	{
	public:
		slimage::ThreadingOptions threadopt;

		Parameters opt;

		ImagePoints points;

		std::vector<Cluster> cluster;

		std::size_t clusterCount() const {
			return cluster.size();
		}

		unsigned int width() const {
			return points.width();
		}

		unsigned int height() const {
			return points.height();
		}

		Clustering();

		std::vector<Seed> getClusterCentersAsSeeds() const;

		void CreatePoints(const slimage::Image3f& image, const slimage::Image1ui16& depth, const slimage::Image3f& normals=slimage::Image3f());

		void CreatePoints(const slimage::Image3ub& image, const slimage::Image1ui16& depth, const slimage::Image3f& normals=slimage::Image3f());

//		/** Find super pixel clusters */
//		void ComputeSuperpixels(const slimage::Image1f& edges);

		std::vector<int> ComputePixelLabels();

		void ComputeSuperpixels(const std::vector<Seed>& seeds);

		slimage::Image1f ComputeDepthDensity();

		std::vector<Seed> FindSeeds();

		std::vector<Seed> FindSeeds(const std::vector<Seed>& old_seeds, const ImagePoints& old_points);

		void ComputeEdges(slimage::Image1f& edges);

		void ImproveSeeds(std::vector<Seed>& seeds, const slimage::Image1f& edges);

		void CreateClusters(const std::vector<Seed>& seeds);

		void MoveClusters();

		SuperpixelGraph CreateNeighborhoodGraph();

		std::vector<unsigned int> CreateSegments(const SuperpixelGraph& neighbourhood, unsigned int* cnt_label);

		/**
		 * Signature of F :
		 * void F(unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p)
		 */
		template<typename F>
		void ForPixelClusters(F f) {
			for(unsigned int i=0; i<cluster.size(); i++) {
				const Cluster& c = cluster[i];
				for(unsigned int p : c.pixel_ids) {
					f(i, c, p, points[p]);
				}
			}
		}

		template<typename F>
		void ForClustersNoReturn(F f) {
			for(unsigned int i=0; i<cluster.size(); i++) {
				f(cluster[i]);
			}
		}

		template<typename F>
		auto ForClusters(F f) -> std::vector<decltype(f(cluster[0]))> {
			std::vector<decltype(f(cluster[0]))> data(cluster.size());
			for(unsigned int i=0; i<cluster.size(); i++) {
				data[i] = f(cluster[i]);
			}
			return data;
		}

		template<typename F>
		auto ForClusterCenters(F f) -> std::vector<decltype(f(cluster[0].center))> {
			std::vector<decltype(f(cluster[0].center))> data(cluster.size());
			for(unsigned int i=0; i<cluster.size(); i++) {
				data[i] = f(cluster[i].center);
			}
			return data;
		}

		void ComputeClusterInfo() {
			return ForClustersNoReturn([this](Cluster& c) { c.ComputeClusterInfo(points, opt); });
		}

		ClusterGroupInfo ComputeClusterGroupInfo();

		inline float DistanceForNormals(const Eigen::Vector3f& x, const Eigen::Vector3f& y)
		{
			// x and y are assumed to be normalized
			// dot(x,y) yields the cos of the angle
			// 1 is perfect, -1 is opposite
			// map to [0|2]
			return 1.0f - x.dot(y);
		}

		template<bool cUseSqrt=false, bool WeightNormalsByDepth=true>
		inline float Distance(const Point& u, const Point& v)
		{
			float d_color;
			if(cUseSqrt) {
				d_color = (u.color - v.color).norm();
				d_color /= 0.25f;
			}
			else {
				d_color = (u.color - v.color).squaredNorm();
				d_color /= (0.25f * 0.25f);
			}
			if(d_color > 1.0f) {
				d_color = 1.0f;
			}

			float d_world;
			if(cUseSqrt) {
				d_world = (u.world - v.world).norm();
				d_world /= opt.base_radius;
			}
			else {
				d_world = (u.world - v.world).squaredNorm();
				d_world /= opt.base_radius * opt.base_radius;
			}

			float d_normal = DistanceForNormals(u.normal, v.normal);
			if(WeightNormalsByDepth) {
				d_normal *= 2.0f / (1.0f + 0.5f * (u.depth() + v.depth()));
			}

			return
				opt.weight_color * d_color
				+ opt.weight_spatial * d_world
				+ opt.weight_normal * d_normal;
		}

		template<bool cUseSqrt=false, bool cDisparity=true, bool WeightNormalsByDepth=true>
		inline float DistancePlanar(const Point& u, const Point& v)
		{
			float d_color;
			if(cUseSqrt) {
				d_color = (u.color - v.color).norm();
				d_color /= 0.25f;
			}
			else {
				d_color = (u.color - v.color).squaredNorm();
				d_color /= (0.25f * 0.25f);
			}
			if(d_color > 1.0f) {
				d_color = 1.0f;
			}

			float d_point;
			if(cUseSqrt) {
				d_point = (u.pos - v.pos).norm();
				d_point /= 0.5f*(u.image_super_radius + v.image_super_radius);
			}
			else {
				d_point = (u.pos - v.pos).squaredNorm();
				float q = 0.5f*(u.image_super_radius + v.image_super_radius);
				d_point /= q*q;
			}

			float d_depth;
			if(cDisparity) {
				if(cUseSqrt) {
					d_depth = std::abs(1.0f/u.depth() - 1.0f/v.depth());
				}
				else {
					d_depth = Square(1.0f/u.depth() - 1.0f/v.depth());
				}
			}
			else {
				if(cUseSqrt) {
					d_depth = std::abs(u.depth() - v.depth());
				}
				else {
					d_depth = Square(u.depth() - v.depth());
				}
			}

			float d_normal = DistanceForNormals(u.normal, v.normal);
			if(WeightNormalsByDepth) {
				d_normal *= 2.0f / (1.0f + 0.5f * (u.depth() + v.depth()));
			}

			return
				opt.weight_color * d_color
				+ opt.weight_spatial * d_point
				+ opt.weight_depth*d_depth
				+ opt.weight_normal * d_normal;
		}

		inline float Distance(const Point& u, const Cluster& c) {
			return Distance(u, c.center);
		}

		inline float Distance(const Cluster& x, const Cluster& y) {
			return Distance(x.center, y.center);
		}

		template<bool cUseSqrt=true>
		inline float DistanceSpatial(const Point& x, const Point& y) {
			float d_point;
			if(cUseSqrt) {
				d_point = (x.pos - y.pos).norm();
			}
			else {
				d_point = (x.pos - y.pos).squaredNorm();
			}
			return d_point;
		}

		inline float DistanceSpatial(const Cluster& x, const Cluster& y) {
			return DistanceSpatial(x.center, y.center);
		}

		template<bool WeightNormalsByDepth=true>
		inline float DistanceWorld(const Cluster& x, const Cluster& y) {
			// mean cluster colors,
			float d_color = (x.center.color - y.center.color).norm();
			// mean cluster normals and
			float d_normal = DistanceForNormals(x.center.normal, y.center.normal);
			if(WeightNormalsByDepth) {
				d_normal *= 2.0f / (1.0f + 0.5f * (x.center.depth() + y.center.depth()));
			}
			// world position of cluster centers.
			float d_world = (x.center.world - y.center.world).norm();
			return opt.weight_color*d_color
				+ opt.weight_spatial*d_world
				+ opt.weight_normal*d_normal;
		}

		template<bool WeightNormalsByDepth=true>
		inline float DistanceColorNormal(const Cluster& x, const Cluster& y) {
			// mean cluster colors,
			float d_color = (x.center.color - y.center.color).norm();
			// mean cluster normals and
			float d_normal = DistanceForNormals(x.center.normal, y.center.normal);
			if(WeightNormalsByDepth) {
				d_normal *= 2.0f / (1.0f + 0.5f * (x.center.depth() + y.center.depth()));
			}
			return opt.weight_color*d_color + opt.weight_normal*d_normal;
		}

		template<bool WeightNormalsByDepth=true>
		inline float DistanceNormal(const Cluster& x, const Cluster& y) {
			// mean cluster normals and
			float d_normal = DistanceForNormals(x.center.normal, y.center.normal);
			if(WeightNormalsByDepth) {
				d_normal *= 2.0f / (1.0f + 0.5f * (x.center.depth() + y.center.depth()));
			}
			return opt.weight_normal*d_normal;
		}

	};

	slimage::Image1f ComputeDepthDensity(const ImagePoints& points, const Parameters& opt);

	slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Seed>& seeds, const slimage::Image1f& target, const Parameters& opt);
}

#endif
