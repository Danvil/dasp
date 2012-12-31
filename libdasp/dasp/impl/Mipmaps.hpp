/*
 * Mipmaps.hpp
 *
 *  Created on: Feb 4, 2012
 *      Author: david
 */

#ifndef SUPERPOINTS_MIPMAPS_HPP_
#define SUPERPOINTS_MIPMAPS_HPP_

#include <Slimage/Slimage.hpp>
#include <Eigen/Dense>
#include <vector>

namespace dasp {
namespace Mipmaps {

Eigen::MatrixXf SumMipMapWithBlackBorder(const Eigen::MatrixXf& img_big);

Eigen::MatrixXf SumMipMap(const Eigen::MatrixXf& img_big);

slimage::Image2f SumMipMapWithAbs(const slimage::Image2f& img_big);

std::vector<Eigen::MatrixXf> ComputeMipmaps(const Eigen::MatrixXf& img, unsigned int min_size);

std::vector<slimage::Image2f> ComputeMipmapsWithAbs(const Eigen::MatrixXf& img, unsigned int min_size);

}}

#endif
