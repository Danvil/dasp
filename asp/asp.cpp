#include "asp.hpp"

namespace asp
{

	namespace impl
	{

		std::vector<std::vector<int>> clusters_from_labels(const std::vector<int>& labels, unsigned int cluster_count)
		{
			std::vector<std::vector<int>> clusters(cluster_count);
			for(unsigned int i=0; i<labels.size(); i++) {
				int label = labels[i];
				if(label >= 0) {
					clusters[label].push_back(i);
				}
			}
			return clusters;
		}

		void remove_empty_clusters(std::vector<std::vector<int>>& clusters)
		{
			auto it = std::remove_if(clusters.begin(), clusters.end(),
				[](const std::vector<int>& c) { return c.empty(); });
			clusters.resize(it - clusters.begin());
		}

		void coordinate_mean(const std::vector<int>& v, int w, float& result_cx, float& result_cy)
		{
			int sx = 0;
			int sy = 0;
			for(int i : v) {
				sx += i % w;
				sy += i / w;
			}
			result_cx = static_cast<float>(sx) / static_cast<float>(v.size());
			result_cy = static_cast<float>(sy) / static_cast<float>(v.size());
		}

	}

	namespace impl
	{
		template<typename F>
		Eigen::Vector3f mean_rgb(const std::vector<F>& f, const std::vector<int>& indices) {
			Eigen::Vector3f sum = Eigen::Vector3f::Zero();
			for(int i : indices) {
				sum += f[i];
			}
			return sum / static_cast<float>(indices.size());
		}
	}

	Parameters parameters_default()
	{
		Parameters p;
		p.num_iterations = 5;
		p.pds_mode = "spds";
		p.seed_mean_radius_factor = 0.25f;
		p.coverage = 3.0f;
		p.weight_compact = 1.0f;
		return p;
	}

	Superpixels<Eigen::Vector3f> AdaptiveSuperpixelsRGB(const Eigen::MatrixXf& density, const std::vector<Eigen::Vector3f>& features, const Parameters& p)
	{
		return AdaptiveSuperpixels(
			// density
			density,
			// feature
			features,
			// metric
			[](const Eigen::Vector3f& a, const Eigen::Vector3f& b) {
				return (a - b).norm();
			},
			// mean
			[](const std::vector<Eigen::Vector3f>& features, const std::vector<int>& indices) {
				return impl::mean_rgb(features, indices);
			},
			// parameters
			p);
	}

}