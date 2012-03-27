/*
 * Superpixels.hpp
 *
 *  Created on: Jan 26, 2012
 *      Author: david
 */

#ifndef SUPERPIXELS_HPP_
#define SUPERPIXELS_HPP_

#include "Tools.hpp"
//#include "SuperpixelHistogram.hpp"
#include "tools/Graph.hpp"
#include <Slimage/Slimage.hpp>
#include <Slimage/Parallel.h>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>

namespace dasp
{
	extern std::map<std::string,slimage::ImagePtr> sDebugImages;

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

		/** circularity = b/a = 1/sqrt(||gradient||^2 + 1) = sqrt(1 - ecc*ecc) */
		float circularity;

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

		/** Computes the normal from the gradient */
		Eigen::Vector3f computeNormal() const {
			return Eigen::Vector3f(gradient[0], gradient[1], -1.0f) * circularity;
		}

//	public:
//		 EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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

	namespace ColorSpaces {
		enum Type {
			RGB, HSV, HN
		};
	}
	typedef ColorSpaces::Type ColorSpace;

	struct Parameters
	{
		Parameters();

		/** camera parameters */
		Camera camera;

		ColorSpace color_space;

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

		/** Desired number of superpixels */
		unsigned int count;

		/** Method used to compute seed points */
		SeedMode seed_mode;

		bool gradient_adaptive_density;

		bool is_conquer_enclaves;

		float segment_threshold;

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
		static constexpr float cPercentage = 0.95f; //0.99f;
		static constexpr float cSigmaScale = 1.959964f; //2.575829f;

		Point center;

		std::vector<unsigned int> pixel_ids;

		bool hasPoints() const {
			return pixel_ids.size() > 3;
		}

		void addPixel(unsigned int index) {
			pixel_ids.push_back(index);
		}

		void removePixel(unsigned int index) {
			auto it = std::find(pixel_ids.begin(), pixel_ids.end(), index);
			if(it != pixel_ids.end()) {
				pixel_ids.erase(it);
			}
		}

//		void addPixels(const std::vector<unsigned int>& v) {
//			pixel_ids.insert(pixel_ids.begin(), v.begin(), v.end());
//		}

		// point covariance matrix
		Eigen::Matrix3f cov;
		// eigenvalues of the covariance matrix
		Eigen::Vector3f ew;
		// eigenvectors of the covariance matrix
		Eigen::Matrix3f ev;

		/** Thickness of the cluster computed using smalles eigenvalue */
		float thickness;
		/** Ratio of plane eigenvalues */
		float circularity;
		/** eccentricity of the ellipse described by a and b */
		float eccentricity;
		/** actual area / expeted area defined by base radius*/
		float area_quotient;

		/** Thickness of cluster computed using orthogonal distance from plane
		 * WARNING: only computed if ComputeExt is called!
		 */
		float thickness_plane;

		/** number of pixel which are within superpixel radius but are not part of the superpixel
		 * WARNING: only computed if ComputeExt is called!
		 */
		float coverage_error;

		/** actual area of the superpixel (computed from all pixels) */
		float area_actual;
		/** expected area of the superpixel (considering local geometry and thickness) */
		float area_expected;
		/** expected area using the actual base radius (computed from cluster count) (same for all clusters...) */
		float area_expected_global;

		void UpdateCenter(const ImagePoints& points, const Parameters& opt);

		void ComputeExt(const ImagePoints& points, const Parameters& opt);

	};

	struct ClusterGroupInfo
	{
		Histogram<float> hist_thickness;
		Histogram<float> hist_circularity;
		Histogram<float> hist_area_quotient;
		Histogram<float> hist_coverage_error;
	};

	void SetRandomNumberSeed(unsigned int seed);

	class Clustering
	{
	public:
		slimage::ThreadingOptions threadopt;

		Parameters opt;

		ImagePoints points;

		slimage::Image1f density;

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

		std::vector<int> ComputePixelLabels() const;

		slimage::Image1i ComputeLabels() const;

		void ComputeSuperpixels(const std::vector<Seed>& seeds);

		void ConquerEnclaves();

		void ConquerMiniEnclaves();

		std::vector<Seed> FindSeeds();

		std::vector<Seed> FindSeeds(const std::vector<Seed>& old_seeds, const ImagePoints& old_points);

		slimage::Image1f ComputeEdges();

		void ImproveSeeds(std::vector<Seed>& seeds, const slimage::Image1f& edges);

		void CreateClusters(const std::vector<Seed>& seeds);

		void PurgeInvalidClusters();

		void MoveClusters();

		std::vector<std::vector<unsigned int> > ComputeBorderPixels(const graph::Graph& graph) const;

		struct NeighborGraphSettings {
			NeighborGraphSettings() {
				cut_by_spatial = true;
				max_spatial_distance_mult = 5.0f;
				min_border_overlap = 0.05f;
				cost_function = NormalColor;
			}
			bool cut_by_spatial;
			float max_spatial_distance_mult;
			float min_border_overlap;
			enum CostFunction {
				SpatialNormalColor,
				NormalColor
			};
			CostFunction cost_function;
		};

		/** Creates the superpixel neighborhood graph. Superpixels are neighbors if they share border pixels. */
		graph::Graph CreateNeighborhoodGraph(NeighborGraphSettings settings=NeighborGraphSettings()) const;

		struct Segmentation
		{
			std::vector<unsigned int> cluster_labels;
			unsigned int segment_count;
			graph::Graph cluster_graph;
			graph::Graph segmentation_graph;

			unsigned int countSegments() const {
				return segment_count;
			}

			unsigned int countClusters() const {
				return cluster_graph.nodes_;
			}
		};

		/** Create clusters of superpixels (= segments) */
		Segmentation CreateSegmentation(const graph::Graph& neighbourhood) const;

		/**
		 * Signature of F :
		 * void F(unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p)
		 */
		template<typename F>
		void ForPixelClusters(F f) const {
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
		void ForClusterCentersNoReturn(F f) {
			for(unsigned int i=0; i<cluster.size(); i++) {
				f(cluster[i].center);
			}
		}

		template<typename F>
		auto ForClusterCenters(F f) -> std::vector<decltype(f(cluster[0].center))> {
			std::vector<decltype(f(cluster[0].center))> data(cluster.size());
			for(unsigned int i=0; i<cluster.size(); i++) {
				data[i] = f(cluster[i].center);
			}
			return data;
		}

		void ComputeExt() {
			return ForClustersNoReturn([this](Cluster& c) { c.ComputeExt(points, opt); });
		}

		ClusterGroupInfo ComputeClusterGroupInfo(unsigned int n, float max_thick);

		template<bool cUseSqrt>
		inline float SpatialDistance(const Point& x, const Point& y) const {
			float d_world;
			if(cUseSqrt) {
				d_world = (x.world - y.world).norm();
				d_world /= opt.base_radius;
			}
			else {
				d_world = (x.world - y.world).squaredNorm();
				d_world /= opt.base_radius * opt.base_radius;
			}
			return d_world;
		}

		template<bool cUseSqrt>
		inline static float ColorDistance(const Point& u, const Point& v) {
			float d_color;
			if(cUseSqrt) {
				d_color = (u.color - v.color).norm();
//				d_color /= 0.25f;
			}
			else {
				d_color = (u.color - v.color).squaredNorm();
//				d_color /= (0.25f * 0.25f);
			}
//			if(d_color > 1.0f) {
//				d_color = 1.0f;
//			}
			return d_color;
		}

		inline static float NormalDistance(const Point& u, const Point& v) {
			// we want to compute 1 - dot(n(u), n(v))
			// the normal is implicitly given by the gradient
			// multiplying with the circularity yields the required normalization
			return 1.0f - (u.gradient.dot(v.gradient) + 1.0f) * u.circularity * v.circularity;
		}

		inline static float DepthWeight(const Point& u, const Point& v) {
			return 2.0f / (1.0f + 0.5f * (u.depth() + v.depth()));
		}

		template<bool cUseSqrt=false, bool WeightNormalsByDepth=false>
		inline float Distance(const Point& u, const Point& v) const
		{
			float d_world = SpatialDistance<cUseSqrt>(u,v);
			float d_color = ColorDistance<cUseSqrt>(u,v);
			float d_normal = NormalDistance(u, v);
			if(WeightNormalsByDepth) {
				d_normal *= DepthWeight(u, v);
			}
			return opt.weight_spatial * d_world
				+ opt.weight_color * d_color
				+ opt.weight_normal * d_normal;
		}

//		template<bool cUseSqrt=true, bool cDisparity=true, bool WeightNormalsByDepth=false>
//		inline float DistancePlanar(const Point& u, const Point& v) const
//		{
//			float d_color = ColorDistance<cUseSqrt>(u,v);
//
//			float d_point;
//			if(cUseSqrt) {
//				d_point = (u.pos - v.pos).norm();
//				d_point /= 0.5f*(u.image_super_radius + v.image_super_radius);
//			}
//			else {
//				d_point = (u.pos - v.pos).squaredNorm();
//				float q = 0.5f*(u.image_super_radius + v.image_super_radius);
//				d_point /= q*q;
//			}
//
//			float d_depth;
//			if(cDisparity) {
//				if(cUseSqrt) {
//					d_depth = std::abs(1.0f/u.depth() - 1.0f/v.depth());
//				}
//				else {
//					d_depth = Square(1.0f/u.depth() - 1.0f/v.depth());
//				}
//			}
//			else {
//				if(cUseSqrt) {
//					d_depth = std::abs(u.depth() - v.depth());
//				}
//				else {
//					d_depth = Square(u.depth() - v.depth());
//				}
//			}
//
//			float d_normal = NormalDistance(u, v);
//			if(WeightNormalsByDepth) {
//				d_normal *= DepthWeight(u, v);
//			}
//
//			return
//				opt.weight_color * d_color
//				+ opt.weight_spatial * d_point
//				+ opt.weight_depth*d_depth
//				+ opt.weight_normal * d_normal;
//		}

//		inline float Distance(const Point& u, const Cluster& c) const {
//			return Distance(u, c.center);
//		}
//
//		inline float Distance(const Cluster& x, const Cluster& y) const {
//			return Distance(x.center, y.center);
//		}

		template<bool WeightNormalsByDepth=false>
		inline float DistanceSpatialColorNormal(const Cluster& x, const Cluster& y) const {
			// world position of cluster centers.
			float d_world = SpatialDistance<true>(x.center, y.center);
			// mean cluster colors,
			float d_color = ColorDistance<true>(x.center, y.center);
			// mean cluster normals and
			float d_normal = NormalDistance(x.center, y.center);
			if(WeightNormalsByDepth) {
				d_normal *= DepthWeight(x.center, y.center);
			}
			return opt.weight_spatial*d_world
				+ opt.weight_color*d_color
				+ opt.weight_normal*d_normal;
		}

		template<bool WeightNormalsByDepth=false>
		inline float DistanceColorNormal(const Cluster& x, const Cluster& y) const {
			// mean cluster colors,
			float d_color = ColorDistance<true>(x.center, y.center);
			// mean cluster normals and
			float d_normal = NormalDistance(x.center, y.center);
			if(WeightNormalsByDepth) {
				d_normal *= DepthWeight(x.center, y.center);
			}
			return opt.weight_color*d_color
				+ opt.weight_normal*d_normal;
		}

		template<bool WeightNormalsByDepth=false>
		inline float DistanceNormal(const Cluster& x, const Cluster& y) const {
			// mean cluster normals and
			float d_normal = NormalDistance(x.center, y.center);
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
