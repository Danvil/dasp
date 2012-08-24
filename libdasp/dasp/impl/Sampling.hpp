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
#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <vector>

namespace dasp
{

	slimage::Image1f ComputeDepthDensity(const ImagePoints& points, const Parameters& opt);

	slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Seed>& seeds, const slimage::Image1f& target);

	slimage::Image1f ComputeDepthDensityFromSeeds(const std::vector<Eigen::Vector2f>& seeds, const slimage::Image1f& target);

	std::vector<Seed> FindSeedsGrid(const ImagePoints& points, const Parameters& opt);

	std::vector<Seed> FindSeedsDepthMipmap(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt);

	std::vector<Seed> FindSeedsDepthBlue(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt);

	std::vector<Seed> FindSeedsDepthFloyd(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt);

	std::vector<Seed> FindSeedsDepthFloydExpo(const ImagePoints& points, const slimage::Image1f& density, const Parameters& opt);

	std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const slimage::Image1f& density_delta, bool delete_small_scala_seeds);

	std::vector<Seed> FindSeedsDelta(const ImagePoints& points, const std::vector<Seed>& old_seeds, const ImagePoints& old_points, const slimage::Image1f& density_new, const Parameters& opt);

}

#endif
