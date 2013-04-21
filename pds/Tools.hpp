#ifndef INCLUDED_PDS_TOOLS_HPP
#define INCLUDED_PDS_TOOLS_HPP

#include <boost/random.hpp>

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
		
		inline Eigen::Vector2f RandomCellPoint(int scale, int x, int y, float gamma)
		{
			float sf = static_cast<float>(scale);
			float xf = static_cast<float>(x);
			float yf = static_cast<float>(y);
			boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > delta(
					Rnd(), boost::uniform_real<float>(0.5f-gamma, 0.5f+gamma));
			return Eigen::Vector2f(sf*(xf + delta()), sf*(yf + delta()));
		}

		inline void ScalePoints(std::vector<Eigen::Vector2f>& pnts, float scale)
		{
			for(Eigen::Vector2f& u : pnts) {
				u *= scale;
			}
		}

	}

}

#endif
