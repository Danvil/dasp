/*
 * Mipmap.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: david
 */

#include "Mipmaps.hpp"
#include <Danvil/Tools/MoreMath.h>
//----------------------------------------------------------------------------//
namespace dasp {
namespace Mipmaps {
//----------------------------------------------------------------------------//

slimage::Image1f SumMipMapWithBlackBorder(const slimage::Image1f& img_big)
{
	size_t w_big = img_big.width();
	size_t h_big = img_big.height();
	// the computed mipmap will have 2^i size
	unsigned int size = Danvil::MoreMath::P2Ceil(std::max(w_big, h_big));
	slimage::Image1f img_small(size / 2, size / 2);
	img_small.fill({0.0f});
	// only the part where at least one of the four pixels lies in the big image is iterated
	// the rest was set to 0 with the fill op
	size_t w_small = w_big / 2 + ((w_big % 2 == 0) ? 0 : 1);
	size_t h_small = h_big / 2 + ((h_big % 2 == 0) ? 0 : 1);
	for(size_t y = 0; y < h_small; y++) {
		size_t y_big = y * 2;
		for(size_t x = 0; x < w_small; x++) {
			size_t x_big = x * 2;
			// We sum over all four pixels in the big image (if they are valid).
			// May by invalid because the big image is considered to be enlarged
			// to have a size of 2^i.
			float sum = 0.0f;
			// Since we only test the part where at least one pixel is in also in the big image
			// we do not need to test that (x_big,y_big) is a valid pixel in the big image.
			const float* p_big = img_big.pointer(x_big, y_big);
			sum += *(p_big);
			if(x_big + 1 < w_big) {
				sum += *(p_big + 1);
			}
			if(y_big + 1 < h_big) {
				sum += *(p_big + w_big);
				if(x_big + 1 < w_big) {
					sum += *(p_big + w_big + 1);
				}
			}
			img_small(x, y) = sum;
		}
	}
	return img_small;
}

slimage::Image1f SumMipMap(const slimage::Image1f& img_big)
{
	size_t w_big = img_big.width();
	size_t h_big = img_big.height();
	// the computed mipmap will have 2^i size
	unsigned int size = Danvil::MoreMath::P2Ceil(std::max(w_big, h_big));
	assert(size == w_big && size == h_big && "SumMipMap: Size must be 2^i!");
	size /= 2;
	slimage::Image1f img_small(size, size);
	for(size_t y = 0; y < size; y++) {
		size_t y_big = y * 2;
		for(size_t x = 0; x < size; x++) {
			size_t x_big = x * 2;
			// We sum over all four corresponding pixels in the big image.
			const float* p_big = img_big.pointer(x_big, y_big);
			float sum = *(p_big) + *(p_big + 1) + *(p_big + h_big) + *(p_big + h_big + 1);
			img_small(x, y) = sum;
		}
	}
	return img_small;
}

slimage::Image2f SumMipMapWithAbs(const slimage::Image2f& img_big)
{
	size_t w_big = img_big.width();
	size_t h_big = img_big.height();
	// the computed mipmap will have 2^i size
	unsigned int size = Danvil::MoreMath::P2Ceil(std::max(w_big, h_big));
	assert(size == w_big && size == h_big && "SumMipMap: Size must be 2^i!");
	size /= 2;
	slimage::Image2f img_small(size, size);
	for(size_t y = 0; y < size; y++) {
		size_t y_big = y * 2;
		for(size_t x = 0; x < size; x++) {
			size_t x_big = x * 2;
			// We sum over all four corresponding pixels in the big image.
			const float* p_big = img_big.pointer(x_big, y_big);
			float sum = *(p_big) + *(p_big + 2) + *(p_big + 2*h_big) + *(p_big + 2*h_big + 2);
			float abs = *(p_big+1) + *(p_big + 3) + *(p_big + 2*h_big+1) + *(p_big + 2*h_big + 3);
			img_small(x, y) = slimage::Pixel2f{sum, abs};
		}
	}
	return img_small;
}

std::vector<slimage::Image1f> ComputeMipmaps(const slimage::Image1f& img, unsigned int min_size)
{
	// find number of required mipmap level
	unsigned int max_size = std::max(img.width(), img.height());
	int n_mipmaps = Danvil::MoreMath::PowerOfTwoExponent(max_size);
	n_mipmaps -= Danvil::MoreMath::PowerOfTwoExponent(min_size);
	BOOST_ASSERT(n_mipmaps >= 1);
	std::vector<slimage::Image1f> mipmaps(n_mipmaps + 1);
	mipmaps[0] = img;
	mipmaps[1] = SumMipMapWithBlackBorder(img);
	// create remaining mipmaps
	for(unsigned int i=2; i<=n_mipmaps; i++) {
		assert(mipmaps[i-1].width() == mipmaps[i-1].height());
		assert(mipmaps[i-1].width() >= 1);
		mipmaps[i] = SumMipMap(mipmaps[i - 1]);
	}
	return mipmaps;
}

slimage::Image2f CreateWithAbs(const slimage::Image1f& img)
{
	slimage::Image2f img_withabs(img.width(), img.height());
	for(unsigned int i=0; i<img.size(); i++) {
		img_withabs(i) = slimage::Pixel2f{img[i], std::abs(img[i])};
	}
	return img_withabs;
}

std::vector<slimage::Image2f> ComputeMipmapsWithAbs(const slimage::Image1f& img, unsigned int min_size)
{
	// find number of required mipmap level
	unsigned int max_size = std::max(img.width(), img.height());
	int n_mipmaps = Danvil::MoreMath::PowerOfTwoExponent(max_size);
	n_mipmaps -= Danvil::MoreMath::PowerOfTwoExponent(min_size);
	BOOST_ASSERT(n_mipmaps >= 1);
	std::vector<slimage::Image2f> mipmaps(n_mipmaps + 1);
	mipmaps[0] = CreateWithAbs(img);
	slimage::Image1f last_mm = SumMipMapWithBlackBorder(img);
	mipmaps[1] = CreateWithAbs(last_mm);
	assert(n_mipmaps >= 4);
	for(unsigned int i=2; i<=4; i++) {
		assert(mipmaps[i-1].width() == mipmaps[i-1].height());
		assert(mipmaps[i-1].width() >= 1);
		last_mm = SumMipMap(last_mm);
		mipmaps[i] = CreateWithAbs(last_mm);
	}
	// create remaining mipmaps
	for(unsigned int i=4; i<=n_mipmaps; i++) {
		assert(mipmaps[i-1].width() == mipmaps[i-1].height());
		assert(mipmaps[i-1].width() >= 1);
		mipmaps[i] = SumMipMapWithAbs(mipmaps[i - 1]);
	}
	return mipmaps;
}

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
