/*
 * Mipmaps.hpp
 *
 *  Created on: Feb 4, 2012
 *      Author: david
 */

#ifndef SUPERPOINTS_MIPMAPS_HPP_
#define SUPERPOINTS_MIPMAPS_HPP_

#include <Slimage/Slimage.hpp>
#include <vector>

namespace dasp {
namespace Mipmaps {

slimage::Image1f SumMipMapWithBlackBorder(const slimage::Image1f& img_big);

slimage::Image1f SumMipMap(const slimage::Image1f& img_big);

std::vector<slimage::Image1f> ComputeMipmaps(const slimage::Image1f& img, unsigned int min_size);

}}

#endif
