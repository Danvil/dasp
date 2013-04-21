#ifndef INCLUDED_PDS_TOOLS_HPP
#define INCLUDED_PDS_TOOLS_HPP

#include <boost/random.hpp>

namespace pds
{

	namespace impl
	{
		inline boost::mt19937& Rnd()
		{
			static boost::mt19937 rnd;
			return rnd;
		}

		inline void RndSeed(unsigned int x)
		{
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

		/** Randomly rounds a float up or down s.t. the expected value is the given value */
		inline unsigned int RandomRound(float x)
		{
			if(x <= 0.0f) {
				return 0;
			}
			float a = std::floor(x);
			float r = x - a;
			boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > uniform01(
					Rnd(), boost::uniform_real<float>(0.0f,1.0f));
			return a + (uniform01() >= r ? 0.0f : 1.0f);
		}

		inline std::vector<unsigned int> RandomSample(const std::vector<float>& v, unsigned int num)
		{
			std::vector<float> a(v.size());
			std::partial_sum(v.begin(), v.end(), a.begin());
//			std::copy(a.begin(), a.end(), std::ostream_iterator<float>(std::cout, ", "));
			float ws = a.back();
			boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > rndv(
					Rnd(), boost::uniform_real<float>(0.0f, ws));
			std::vector<unsigned int> idx(num);
			std::generate(idx.begin(), idx.end(),
				[&rndv,&a]() -> unsigned int {
					float x = rndv();
					auto it = std::lower_bound(a.begin(), a.end(), x);
					if(it == a.end()) {
						return a.size() - 1;
					}
					else {
						return it - a.begin();
					}
				});
			return idx;
		}

	}

}

#endif
