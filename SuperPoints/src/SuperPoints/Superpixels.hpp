/*
 * Superpixels.hpp
 *
 *  Created on: Jan 26, 2012
 *      Author: david
 */

#ifndef SUPERPIXELS_HPP_
#define SUPERPIXELS_HPP_

#include <Slimage/Slimage.hpp>
#include <Slimage/Parallel.h>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>

/** Depth-adaptive super pixels */
namespace dasp
{
	typedef float S;

	struct Position
	{
		int x, y;
	};

	struct Seed
	{
		int x, y;
		float scala;
	};

	struct Point
	{
//			bool valid_;
		Eigen::Vector3f color;
		Eigen::Vector2f pos;
//			Eigen::Vector3f world;
		Eigen::Vector3f normal;

		/// depth of point
		S depth;

		/// project size of a super pixel at point depth
		S scala;

		/** Radius of a super pixel at point depth */
		int radius() const {
			return int(0.603553391f * std::round(scala));
			// 0.603553391 = (0.5 + sqrt(2)/2)/2 // TODO what factor?
			// 0.5 would be inner circle, sqrt(2)/2 would be outer circle
		}

		/** Estimated number of super pixels at this point */
		float estimatedCount() const {
			if(scala == 0.0f) {
				return 0.0f;
			}
			// scala * scala is size of a projected super pixel at this depth.
			// so the inverse is the estimated number of super pixels
			return 1.0f / (scala * scala);
		}

		bool valid() const {
			return depth > 0.01f;
		}

		int spatial_x() const {
			return int(pos[0]);
		}

		int spatial_y() const {
			return int(pos[1]);
		}

		static Point Zero() {
			Point p;
//				p.valid_ = false;
			p.color = Eigen::Vector3f::Zero();
			p.pos = Eigen::Vector2f::Zero();
//				p.world = Eigen::Vector3f::Zero();
			p.normal = Eigen::Vector3f::Zero();
			p.depth = 0.0f;
//				p.scala = 0.0f;
			return p;
		}

		Point& operator+=(const Point& x) {
			color += x.color;
			pos += x.pos;
//				world += x.world;
			normal += x.normal;
			depth += x.depth;
//				scala += x.scala;
			return *this;
		}

		Point& operator*=(S s) {
			color *= s;
			pos *= s;
//				world *= s;
			normal *= s;
			depth *= s;
//				scala *= s;
			return *this;
		}

	public:
		 EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	};

	namespace SeedModes {
		enum Type {
			EquiDistant,
			DepthDependentShooting,
			DepthDependentMipmap,
			BlueNoise
		};
	}
	typedef SeedModes::Type SeedMode;

	struct Parameters
	{
		/// number of clusters
		unsigned int cluster_count;

		S weight_spatial;
//			S weight_world;
		S weight_normal;
		S weight_depth;

		/// number of iterations
		unsigned int iterations;

		/// scale for cluster scala for search area
		float coverage;

		/// camera focal length
		float focal;

		/// superpixels per m^2
		float roh;

		SeedMode seed_mode;

		/**
		 * n = roh * (z/f)^2
		 */
		float computeClusterCountFactor() const {
			return roh / (focal * focal);
		}

		/**
		 * s_p = s_r * roh / z = f / sqrt(roh) / z
		 */
		float computePixelSizeFactor() const {
			return focal / std::sqrt(roh);
		}
	};

	struct ParametersExt
	: public Parameters
	  {
		ParametersExt() {}
		ParametersExt(const Parameters& p) : Parameters(p) {}
		unsigned int width, height;
		unsigned int cluster_nx, cluster_ny;
		unsigned int cluster_dx, cluster_dy;
		float radius;
		float spatial_normalizer;
		float weight_spatial_final;
	};

	template<typename K>
	K Square(K x) { return x*x; }

//		inline S Distance(const Point& u, const Point& v, const ParametersExt& opt) {
//			S d_color = (u.color - v.color).squaredNorm();
//			S d_point = S(Square(u.spatial.x - v.spatial.x)) + S(Square(u.spatial.y - v.spatial.y));
//			S d_depth = Square(u.depth - v.depth);
//			return d_color + opt.weight_spatial_final*d_point + opt.weight_depth*d_depth;
//			//return std::sqrt(d_color) + opt.weight_spatial*std::sqrt(d_point) + opt.weight_depth*std::sqrt(d_depth); // TODO more accurate
//		}

	struct ImagePoints {
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

	struct Cluster
	{
		Point center;

		std::vector<unsigned int> pixel_ids;

		bool is_valid() const {
			return pixel_ids.size() > 0;
		}

		void UpdateCenter(const ImagePoints& points);

	};

	inline S DistanceForNormals(const Eigen::Vector3f& x, const Eigen::Vector3f& y) {
		// x and y are assumed to be normalized
		// dot(x,y) yields the cos of the angle
		// 1 is perfect, -1 is opposite
		// map to [0|2]
		return 1.0f - x.dot(y);
		//return -0.434783f + 1.0f/(1.3f + x.dot(y));
	}

	template<bool cUseSqrt=true, bool cDisparity=true>
	inline S Distance(const Point& u, const Point& v, const ParametersExt& opt) {
		S d_color, d_point;
		if(cUseSqrt) {
			d_color = (u.color - v.color).norm();
			d_point = (u.pos - v.pos).norm();
//				d_world = (u.world - v.world).norm();
		}
		else {
			d_color = (u.color - v.color).squaredNorm();
			d_point = (u.pos - v.pos).squaredNorm();
//				d_world = (u.world - v.world).squaredNorm();
		}

//			float mean_depth = 0.5f*(u.depth + v.depth);
//			d_world /= mean_depth * mean_depth;
//			d_point /= mean_depth;

		S d_depth;
//			if(u.valid() && v.valid()) {
		if(cDisparity) {
			if(cUseSqrt) {
				d_depth = std::abs(1.0f/u.depth - 1.0f/v.depth);
			}
			else {
				d_depth = Square(1.0f/u.depth - 1.0f/v.depth);
			}
		}
		else {
			if(cUseSqrt) {
				d_depth = std::abs(u.depth - v.depth);
			}
			else {
				d_depth = Square(u.depth - v.depth);
			}
		}
//			}
//			else {
//				// use some kind of max distance for world, normal and depth values
//				d_depth = 10.0f;
//			}
		S d_normal = DistanceForNormals(u.normal, v.normal);
		// both points have valid depth information
		return d_color
			+ opt.weight_spatial_final*d_point
//				+ opt.weight_world*d_world
			+ opt.weight_normal*d_normal
			+ opt.weight_depth*d_depth;
	}

	inline S Distance(const Point& u, const Cluster& c, const ParametersExt& opt) {
		return Distance(u, c.center, opt);
	}

	inline S Distance(const Cluster& x, const Cluster& y, const ParametersExt& opt) {
		return Distance(x.center, y.center, opt);
	}

	template<bool cUseSqrt=true>
	inline S DistanceSpatial(const Point& x, const Point& y) {
		S d_point;
		if(cUseSqrt) {
			d_point = (x.pos - y.pos).norm();
		}
		else {
			d_point = (x.pos - y.pos).squaredNorm();
		}
		return d_point;
	}

	inline S DistanceSpatial(const Cluster& x, const Cluster& y) {
		return DistanceSpatial(x.center, y.center);
	}

	inline S DistanceWorld(const Cluster& x, const Cluster& y, const ParametersExt& opt) {
		// TODO compute distance between:
		// mean cluster colors,
		S d_color = (x.center.color - y.center.color).norm();
		// mean cluster normals and
		S d_normal = DistanceForNormals(x.center.normal, y.center.normal);
		// world position of cluster centers.
		Eigen::Vector3f pw1(
				x.center.depth*(float(x.center.pos[0]) - float(640/2))/580.0f,
				x.center.depth*(float(x.center.pos[1]) - float(480/2))/580.0f,
				x.center.depth);
		Eigen::Vector3f pw2(
				y.center.depth*(float(y.center.pos[0]) - float(640/2))/580.0f,
				y.center.depth*(float(y.center.pos[1]) - float(480/2))/580.0f,
				y.center.depth);
		S d_world = (pw1 - pw2).norm();
		return d_color
			+ opt.weight_spatial*d_world
			+ opt.weight_normal*d_normal;
		return DistanceSpatial(x.center, y.center);
	}

	ParametersExt ComputeParameters(const Parameters& opt, unsigned int width, unsigned int height);

	ImagePoints CreatePoints(
			const slimage::Image3ub& image,
			const slimage::Image1ui16& depth,
//				const Danvil::Images::Image3fPtr& points,
			const slimage::Image3f& normals,
			const ParametersExt& opt
			);

	/** Find super pixel clusters */
	std::vector<Cluster> ComputeSuperpixels(const ImagePoints& points, const slimage::Image1f& edges, const ParametersExt& opt);

	std::vector<Cluster> ComputeSuperpixels(const ImagePoints& points, const std::vector<Seed>& seeds, const ParametersExt& opt);

	std::vector<Seed> FindSeeds(const ImagePoints& points, const ParametersExt& opt);

	void ComputeEdges(const ImagePoints& points, slimage::Image1f& edges, const ParametersExt& opt,slimage::ThreadingOptions threadopt);

	void ImproveSeeds(std::vector<Seed>& seeds, const ImagePoints& points, const slimage::Image1f& edges, const ParametersExt& opt);

	std::vector<Cluster> CreateClusters(const std::vector<Seed>& seeds, const ImagePoints& points, const ParametersExt& opt);

	void MoveClusters(std::vector<Cluster>& clusters, const ImagePoints& points, const ParametersExt& opt);

	void PlotCluster(const Cluster& cluster, const ImagePoints& points, const slimage::Image3ub& img);

	void PlotCluster(const std::vector<Cluster>& clusters, const ImagePoints& points, const slimage::Image3ub& img);

	void PlotEdges(const std::vector<int>& point_labels, const slimage::Image3ub& img, unsigned int edge_w, unsigned char edge_r, unsigned char edge_g, unsigned char edge_b);

	template<typename F>
	void ForPixelClusters(const std::vector<Cluster>& clusters, const ImagePoints& points, F f) {
		for(unsigned int i=0; i<clusters.size(); i++) {
			const Cluster& cluster = clusters[i];
			for(unsigned int p : cluster.pixel_ids) {
				f(i, cluster, p, points[p]);
			}
		}
	}

	template<typename F>
	std::vector<float> ClassifyClusters(const std::vector<Cluster>& clusters, F f) {
		std::vector<float> values(clusters.size());
		for(unsigned int i=0; i<clusters.size(); i++) {
			values[i] = f(clusters[i].center);
		}
		return values;
	}

}

#endif
