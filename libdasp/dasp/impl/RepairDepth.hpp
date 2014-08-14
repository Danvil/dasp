/*
 * RepairDepth.hpp
 *
 *  Created on: Mar 27, 2012
 *      Author: david
 */

#ifndef DASP_REPAIRDEPTH_HPP_
#define DASP_REPAIRDEPTH_HPP_

#include <slimage/image.hpp>

namespace dasp
{

	void RepairDepth(slimage::Image1ui16& depth, const slimage::Image3ub& color);

	void SmoothDepth(slimage::Image1ui16& depth, const slimage::Image3ub& color);

}

#endif
