/*
 * Sampling.hpp
 *
 *  Created on: Aug 24, 2012
 *      Author: david
 */

#ifndef DASP_IMPL_SAMPLING_HPP_
#define DASP_IMPL_SAMPLING_HPP_

#include "../Point.hpp"
#include "../Seed.hpp"
#include <Eigen/Dense>
#include <vector>

namespace dasp
{
	Eigen::MatrixXf ComputeDepthDensity(const ImagePoints& points, const Parameters& opt);

	Eigen::MatrixXf ComputeDepthDensityFromSeeds(const std::vector<Seed>& seeds, const Eigen::MatrixXf& target);

	Eigen::MatrixXf ComputeDepthDensityFromSeeds(const std::vector<Eigen::Vector2f>& seeds, const Eigen::MatrixXf& target);

	Eigen::MatrixXf ComputeSaliency(const ImagePoints& points, const Parameters& opt);

	void AdaptClusterRadiusBySaliency(ImagePoints& points, const Eigen::MatrixXf& saliency, const Parameters& opt);

	std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const Eigen::MatrixXf& density_delta, bool delete_small_scala_seeds);

	std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const ImagePoints& old_points, const Eigen::MatrixXf& density_new);

}

#endif
