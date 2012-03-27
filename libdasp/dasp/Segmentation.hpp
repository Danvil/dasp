/*
 * Segmentation.hpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

#ifndef SEGMENTATION_HPP_
#define SEGMENTATION_HPP_

#include <Slimage/Slimage.hpp>
#include "Superpixels.hpp"
#include <vector>

namespace dasp
{

extern std::vector<slimage::Image3ub> cSegmentationDebug;

slimage::Image1ub ComputeBoundary(const Clustering& clusters);

}

#endif
