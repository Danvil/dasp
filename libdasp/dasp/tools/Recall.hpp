/*
 * Recall.hpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#ifndef DASP_RECALL_HPP_
#define DASP_RECALL_HPP_

#include <Slimage/Slimage.hpp>

namespace dasp
{

float ComputeRecallBox(const slimage::Image1ub& img_exp, const slimage::Image1ub& img_act, int d);

float ComputeRecallGaussian(const slimage::Image1ub& img_exp, const slimage::Image1ub& img_act, float sigma);

}

#endif
