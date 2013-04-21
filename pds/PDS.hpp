#ifndef INCLUDED_PDS_PDS_HPP
#define INCLUDED_PDS_PDS_HPP

#include <Eigen/Dense>
#include <boost/random.hpp>
#include <vector>

namespace pds
{

	namespace impl
	{
		inline boost::mt19937& Rnd() {
			static boost::mt19937 rnd;
			return rnd;
		}

		inline void RndSeed(unsigned int x) {
			Rnd().seed(x);
		}
	}

	std::vector<Eigen::Vector2f> RectGrid(const Eigen::MatrixXf& density);

	std::vector<Eigen::Vector2f> HexGrid(const Eigen::MatrixXf& density);

	std::vector<Eigen::Vector2f> SimplifiedPoissonDiscSamplingOld(const Eigen::MatrixXf& density);

	std::vector<Eigen::Vector2f> SimplifiedPoissonDiscSampling(const Eigen::MatrixXf& density);

	std::vector<Eigen::Vector2f> Fattal(const Eigen::MatrixXf& density, unsigned int max_steps=0);

}

#endif
