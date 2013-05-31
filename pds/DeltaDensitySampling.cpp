
#include "PDS.hpp"
#include "Tools.hpp"
#include <density/ScalePyramid.hpp>
#include <density/PointDensity.hpp>
#include <iostream>

#define VERBOSE

namespace pds
{

	std::vector<Eigen::Vector2f> DeltaDensitySampling(const Eigen::MatrixXf& dnew, const std::vector<Eigen::Vector2f>& seeds_old)
	{
	#ifdef CREATE_DEBUG_IMAGES
		slimage::Image3ub debug(points.width(), points.height(), {{0,0,0}});
		sDebugImages["seeds_delta"] = slimage::Ptr(debug);
	#endif
		// compute old density
		Eigen::MatrixXf dprev = density::PointDensity(seeds_old, dnew);
#ifdef VERBOSE
		std::cout << "DDS: seeds_old.size()=" << seeds_old.size() << std::endl;
#endif
		// delta density
		Eigen::MatrixXf ddelta = dnew - dprev;
		Eigen::MatrixXf ddelta_abs = ddelta.cwiseAbs();
		Eigen::MatrixXf dsub = 0.5f * (ddelta_abs - ddelta);
#ifdef VERBOSE
		std::cout << "DDS: dsub.sum()=" << dsub.sum() << std::endl;
#endif
		Eigen::MatrixXf dadd = 0.5f * (ddelta_abs + ddelta);
#ifdef VERBOSE
		std::cout << "DDS: dadd.sum()=" << dadd.sum() << std::endl;
#endif
		// sample points
		std::vector<Eigen::Vector2f> psub = SimplifiedPDS(dsub);
#ifdef VERBOSE
		std::cout << "DDS: psub.size()=" << psub.size() << std::endl;
#endif
		std::vector<Eigen::Vector2f> padd = SimplifiedPDS(dadd);
#ifdef VERBOSE
		std::cout << "DDS: padd.size()=" << padd.size() << std::endl;
#endif
		// remove points
		// TODO: this is O(N*M) and slow
		std::vector<Eigen::Vector2f> samples = seeds_old;
		for(const Eigen::Vector2f& p : psub) {
			float d_min = (p - samples.front()).squaredNorm();
			auto it_min = samples.begin();
			for(auto it=samples.begin(); it!=samples.end(); ++it) {
				float d = (p - *it).squaredNorm();
				if(d < d_min) {
					d_min = d;
					it_min = it;
				}
			}
			samples.erase(it_min);
		}
		// add points
		samples.insert(samples.end(), padd.begin(), padd.end());
#ifdef VERBOSE
		std::cout << "DDS: samples.size()=" << samples.size() << std::endl;
#endif
		return samples;
	}

}
