/*
 * Recall.hpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#include "Recall.hpp"
#include <cmath>
#include <boost/assert.hpp>

namespace dasp
{

const unsigned int cBorder = 50;

template<typename K, typename L>
float ComputeRecallBoxImpl(const slimage::Image<slimage::Traits<K,1>>& img_exp, const slimage::Image<slimage::Traits<L,1>>& img_act, K threshold_exp, L threshold_act, int d)
{
	BOOST_ASSERT(img_exp.hasSameShape(img_act));
	// check how much pixels from the expected boundary are near a boundary pixel in the actual image
	unsigned int cnt = 0;
	unsigned int cnt_recalled = 0;
	for(int y=cBorder+d; y+cBorder+d<static_cast<int>(img_exp.height()); y++) {
		for(int x=cBorder+d; x+cBorder+d<static_cast<int>(img_exp.width()); x++) {
			if(img_exp(x,y) >= threshold_act) {
				bool recalled = false;
				for(int u=-d; u<=+d; u++) {
					for(int v=-d; v<=+d; v++) {
						if(img_act(x+u, y+v) >= threshold_exp) {
							recalled = true;
						}
					}
				}
				cnt ++;
				if(recalled) {
					cnt_recalled ++;
				}
			}
		}
	}
	return cnt == 0 ? 1.0f : static_cast<float>(cnt_recalled) / static_cast<float>(cnt);
}

float ComputeRecallBox(const slimage::Image1ub& img_exp, const slimage::Image1ub& img_act, int d)
{
	return ComputeRecallBoxImpl(img_exp, img_act, static_cast<unsigned char>(255), static_cast<unsigned char>(255), d);
}

PrecisionRecallStats PrecisionRecall(const slimage::Image1ub& img_relevant, const slimage::Image1ub& img_retrievend, int d)
{
	BOOST_ASSERT(img_relevant.hasSameShape(img_retrievend));
	PrecisionRecallStats stats;
	for(int y=cBorder+d; y+cBorder+d<static_cast<int>(img_relevant.height()); y++) {
		for(int x=cBorder+d; x+cBorder+d<static_cast<int>(img_relevant.width()); x++) {
			if(img_relevant(x,y)) {
				stats.num_relevant ++;
				bool is_retrieved = false;
				// check if pixels from the relevant boundary is near a boundary pixel in the retrieved image
				for(int u=-d; u<=+d && !is_retrieved; u++) {
					for(int v=-d; v<=+d && !is_retrieved; v++) {
						if(img_retrievend(x+u, y+v)) {
							is_retrieved = true;
						}
					}
				}
				if(is_retrieved) {
					stats.num_relevant_and_retrieved ++;
				}
			}
			if(img_retrievend(x,y)) {
				stats.num_retrieved ++;
				bool is_relevant = false;
				// check if pixels from the retrieved boundary is near a boundary pixel in the relevant image
				for(int u=-d; u<=+d && !is_relevant; u++) {
					for(int v=-d; v<=+d && !is_relevant; v++) {
						if(img_relevant(x+u, y+v)) {
							is_relevant = true;
						}
					}
				}
				if(is_relevant) {
					stats.num_retrieved_and_relevant ++;
				}
			}
		}
	}
	return stats;
}

float ComputeRecallGaussian(const slimage::Image1ub& img_exp, const slimage::Image1ub& img_act, float sigma)
{
	BOOST_ASSERT(img_exp.hasSameShape(img_act));
	const float exp_arg_norm = -0.5f / (sigma * sigma);
	const float cResponseThreshold = 0.05f;
	int d = std::round(sigma * std::sqrt( - 2.0f * std::log(cResponseThreshold)));
	// check how much pixels from the expected boundary are near a boundary pixel in the actual image
	unsigned int cnt = 0;
	float recalled = 0;
	for(int y=cBorder+d; y+cBorder+d<static_cast<int>(img_exp.height()); y++) {
		for(int x=cBorder+d; x+cBorder+d<static_cast<int>(img_exp.width()); x++) {
			if(img_exp(x,y)) {
				float d2_min = 1e9f;
				for(int u=-d; u<=+d; u++) {
					for(int v=-d; v<=+d; v++) {
						if(img_act(x+u, y+v)) {
							float d2 = static_cast<float>(u*u + v*v);
							if(d2 < d2_min) {
								d2_min = d2;
							}
						}
					}
				}
				cnt ++;
				recalled += std::exp(exp_arg_norm*d2_min);
			}
		}
	}
	return recalled / static_cast<float>(cnt);
}

}
