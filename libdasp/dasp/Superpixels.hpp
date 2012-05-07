/*
 * Superpixels.hpp
 *
 *  Created on: Jan 26, 2012
 *      Author: david
 */

#ifndef SUPERPIXELS_HPP_
#define SUPERPIXELS_HPP_

#include "Point.hpp"
#include "Clustering.hpp"
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

	struct ClusterGroupInfo
	{
		Histogram<float> hist_thickness;
		Histogram<float> hist_circularity;
		Histogram<float> hist_area_quotient;
		Histogram<float> hist_coverage_error;
	};

	void SetRandomNumberSeed(unsigned int seed);

	class Superpixels
	{
	public:
		slimage::ThreadingOptions threadopt;

		Parameters opt;

		slimage::Image3ub color_raw;

		ImagePoints points;

		slimage::Image1f density;

		std::vector<Cluster> cluster;

		std::vector<Seed> seeds_previous;
		std::vector<Seed> seeds;

		std::size_t clusterCount() const {
			return cluster.size();
		}

		unsigned int width() const {
			return points.width();
		}

		unsigned int height() const {
			return points.height();
		}

		Superpixels();

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

		std::vector<Seed> FindSeeds(const ImagePoints& old_points);

		slimage::Image1f ComputeEdges();

		void ImproveSeeds(std::vector<Seed>& seeds, const slimage::Image1f& edges);

		void CreateClusters(const std::vector<Seed>& seeds);

		void PurgeInvalidClusters();

		void MoveClusters();

		std::vector<std::vector<unsigned int> > ComputeBorderPixels(const graph::Graph& graph) const;

		/** Computes points which lie on the border between segments
		 * @param list of point indices
		 */
		std::vector<unsigned int> ComputeBorderPixelsComplete() const;

		struct NeighborGraphSettings {
			NeighborGraphSettings() {
				cut_by_spatial = true;
				max_spatial_distance_mult = 5.0f;
				min_border_overlap = 0.00f;
				min_abs_border_overlap = 1;
				cost_function = NormalColor;
			}
			bool cut_by_spatial;
			float max_spatial_distance_mult;
			float min_border_overlap;
			unsigned min_abs_border_overlap;
			enum CostFunction {
				SpatialNormalColor,
				NormalColor,
				Color
			};
			CostFunction cost_function;
		};

		/** Creates the superpixel neighborhood graph. Superpixels are neighbors if they share border pixels. */
		graph::Graph CreateNeighborhoodGraph(NeighborGraphSettings settings=NeighborGraphSettings()) const;

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

		Eigen::Vector3f ColorToRGB(const Eigen::Vector3f& c) const;

		Eigen::Vector3f ColorFromRGB(const Eigen::Vector3f& c) const;


	};

	std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const slimage::Image1f& density_delta);

	slimage::Image1f ComputeDepthDensity(const ImagePoints& points, const Parameters& opt);

	slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Seed>& seeds, const slimage::Image1f& target);

	slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Eigen::Vector2f>& seeds, const slimage::Image1f& target);

	Superpixels ComputeSuperpixels(const slimage::Image3ub& color, const slimage::Image1ui16& depth, const Parameters& opt);

	void ComputeSuperpixelsIncremental(Superpixels& clustering, const slimage::Image3ub& color, const slimage::Image1ui16& depth);


}

#endif
