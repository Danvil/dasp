#ifndef INCLUDED_ASP_ASP_HPP
#define INCLUDED_ASP_ASP_HPP

#include <pds/PDS.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>

namespace asp
{

	template<typename F>
	struct Cluster
	{
		float x, y;
		F f;
	};

	namespace impl
	{
		inline float superpixel_radius_px(float density)
		{
			constexpr float PI = 3.1415f;
			if(density == 0.0f) {
				return 0.0f;
			}
			else {
				return 1.0f / std::sqrt(PI * density);
			}
		}

		inline float superpixel_radius_px(const Eigen::MatrixXf& density, float px, float py)
		{
			return superpixel_radius_px(density(
					static_cast<int>(px + 0.5f),
					static_cast<int>(py + 0.5f)));
		}

		template<typename F>
		void iterate_box(int width, int height, int cx, int cy, int r, F f)
		{
			const int xmin = std::max<int>(static_cast<int>(cx - r), 0);
			const int xmax = std::min<int>(static_cast<int>(cx + r), width-1);
			const int ymin = std::max<int>(static_cast<int>(cy - r), 0);
			const int ymax = std::min<int>(static_cast<int>(cy + r), height-1);
			for(int y=ymin; y<=ymax; y++) {
				const int pnt_index_0 = y*width;
				for(int x=xmin; x<=xmax; x++) {
					f(x, y, pnt_index_0 + x);
				}
			}
		} 

		/** Density-Adaptive Local Iterative Clustering (DALIC)
		 * ASSUME: density(x,y) == 0 <=> (x,y) invalid
		 */
		template<typename F, typename METRIC>
		std::vector<int> cluster_pixels(const std::vector<Cluster<F>>& clusters, const Eigen::MatrixXf& density, const std::vector<F>& features, METRIC metric, float p_coverage, float p_weight_compact)
		{
			const int width = density.rows();
			const int height = density.cols();
			std::vector<int> labels(width*height, -1);
			std::vector<float> v_dist(labels.size(), 1e9);
			// for each cluster check possible points
			for(int j=0; j<clusters.size(); j++) {
				const Cluster<F>& c = clusters[j];
				const float cluster_radius_px = superpixel_radius_px(density, c.x, c.y);
				if(cluster_radius_px == 0.0f) {
					continue;
				}
				const float alpha = p_weight_compact / cluster_radius_px;
				const float R = cluster_radius_px * p_coverage;

				// iterate_box(width, height, c.x, c.y, R,
				// 	[j,&c,&features,&v_dist,&labels,alpha](int x, int y, int pnt_index) {
				// 		const float dy = static_cast<float>(y) - c.y;
				// 		const float dx = static_cast<float>(x) - c.x;
				// 		const float dist = alpha*std::sqrt(dx*dx + dy*dy) + metric(features[pnt_index], c.f);
				// 		float& v_dist_best = v_dist[pnt_index];
				// 		if(dist < v_dist_best) {
				// 			v_dist_best = dist;
				// 			labels[pnt_index] = j;
				// 		}
				// 	});

				const int xmin = std::max<int>(static_cast<int>(c.x - R), 0);
				const int xmax = std::min<int>(static_cast<int>(c.x + R), width-1);
				const int ymin = std::max<int>(static_cast<int>(c.y - R), 0);
				const int ymax = std::min<int>(static_cast<int>(c.y + R), height-1);
				for(int y=ymin; y<=ymax; y++) {
					const int pnt_index_0 = y*width;
					const float dy = static_cast<float>(y) - c.y;
					const float dy2 = dy*dy;
					for(int x=xmin; x<=xmax; x++) {
						const int pnt_index = pnt_index_0 + x;
						if(density(x,y) == 0.0f) {
							continue;
						}
						const float dx = static_cast<float>(x) - c.x;
						const float dist = alpha*std::sqrt(dx*dx + dy2) + metric(features[pnt_index], c.f);
						float& v_dist_best = v_dist[pnt_index];
						if(dist < v_dist_best) {
							v_dist_best = dist;
							labels[pnt_index] = j;
						}
					}
				}
			}
			return labels;
		}

		std::vector<std::vector<int>> clusters_from_labels(const std::vector<int>& v, unsigned int cluster_count);

		void remove_empty_clusters(std::vector<std::vector<int>>& clusters);

		void coordinate_mean(const std::vector<int>& v, int w, float& result_cx, float& result_cy);

		template<typename F, typename MEAN>
		F local_mean(const Eigen::MatrixXf& density, const std::vector<F>& features, float px, float py, MEAN mean, float p_scale)
		{
			const int width = density.rows();
			const int height = density.cols();
			const int cx = static_cast<int>(px + 0.5f);
			const int cy = static_cast<int>(py + 0.5f);
			const float cluster_radius_px = impl::superpixel_radius_px(density(cx,cy));
			const int r = static_cast<int>(p_scale * cluster_radius_px);
			std::vector<int> indices;
			indices.reserve(4*(r+1)*(r+1)); // >= (2*r+1)^2
			iterate_box(width, height, cx, cy, r,
				[&indices](int x, int y, int i) {
					indices.push_back(i);
				});
			return mean(features, indices);
		}

	}

	template<typename F>
	struct Superpixels
	{
		std::vector<Cluster<F>> clusters;
		std::vector<std::vector<int>> cluster_indices;
		std::vector<int> labels;
	};

	struct Parameters
	{
		int num_iterations;
		std::string pds_mode;
		float seed_mean_radius_factor;
		float coverage;
		float weight_compact;
	};

	Parameters parameters_default();

	/** Computes Adaptive Superpixels
	 * Syntax:
	 *		FEATURE: int x int -> F
	 *		METRIC: F x F -> real
	 *		MEAN: FEATURE x vector<int> -> F
	 *
	 * @param density density function
	 * @param feature feature vector lookup
	 * @param metric feature vector metric
	 * @param mean feature vector mean
	 * @return pixel cluster index
	 */
	template<typename F, typename METRIC, typename MEAN>
	Superpixels<F> AdaptiveSuperpixels(const Eigen::MatrixXf& density, const std::vector<Eigen::Vector2f>& seeds, const std::vector<F>& features, METRIC metric, MEAN mean, const Parameters& p)
	{
		const int width = density.rows();
		Superpixels<F> sp;
		// seed clusters
		std::vector<Eigen::Vector2f> pnts = pds::PoissonDiscSampling(p.pds_mode, density, seeds);
		sp.clusters.resize(pnts.size());
		std::transform(pnts.begin(), pnts.end(), sp.clusters.begin(),
			[&density, &features, mean, &p](const Eigen::Vector2f& pnt) {
				float px = pnt[0];
				float py = pnt[1];
				return Cluster<F> {
					px, py,
					impl::local_mean(density, features, px, py, mean, p.seed_mean_radius_factor)
				};
			});
		// dalic
		for(int k=0; k<p.num_iterations; k++) {
			// update pixel-cluster association
			sp.labels = impl::cluster_pixels(sp.clusters, density, features, metric, p.coverage, p.weight_compact);
			// update clusters
			sp.cluster_indices = impl::clusters_from_labels(sp.labels, sp.clusters.size());
			impl::remove_empty_clusters(sp.cluster_indices);
			sp.clusters.resize(sp.cluster_indices.size());
			for(std::size_t i=0; i<sp.clusters.size(); i++) {
				auto& c = sp.clusters[i];
				c.f = mean(features, sp.cluster_indices[i]);
				impl::coordinate_mean(sp.cluster_indices[i], width, c.x, c.y);
			}
		}
		// write back new indices
		std::fill(sp.labels.begin(), sp.labels.end(), -1);
		for(int i=0; i<sp.cluster_indices.size(); i++) {
			for(int j : sp.cluster_indices[i]) {
				sp.labels[j] = i;
			}
		}
		return sp;
	}

	Superpixels<Eigen::Vector3f> AdaptiveSuperpixelsRGB(const Eigen::MatrixXf& density, const std::vector<Eigen::Vector2f>& seeds, const std::vector<Eigen::Vector3f>& features, const Parameters& p);

}

#endif
