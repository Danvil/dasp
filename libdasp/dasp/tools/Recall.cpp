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

float ComputeRecallBox(const slimage::Image1ub& img_exp, const slimage::Image1ub& img_act, int d)
{
	BOOST_ASSERT(img_exp.hasSameShape(img_act));
	// check how much pixels from the expected boundary are near a boundary pixel in the actual image
	unsigned int cnt = 0;
	unsigned int cnt_recalled = 0;
	for(int y=cBorder+d; y+cBorder+d<static_cast<int>(img_exp.height()); y++) {
		for(int x=cBorder+d; x+cBorder+d<static_cast<int>(img_exp.width()); x++) {
			if(img_exp(x,y)) {
				bool recalled = false;
				for(int u=-d; u<=+d; u++) {
					for(int v=-d; v<=+d; v++) {
						if(img_act(x+u, y+v)) {
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
	return static_cast<float>(cnt_recalled) / static_cast<float>(cnt);
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
