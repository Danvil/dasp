/*
 * Parameters.hpp
 *
 *  Created on: Apr 4, 2012
 *      Author: david
 */

#ifndef DASP_PARAMETERS_HPP_
#define DASP_PARAMETERS_HPP_

#include "Tools.hpp"
#include <ctype.h>

namespace dasp
{

	namespace SeedModes
	{
		enum Type {
			EquiDistant,
			DepthShooting,
			DepthMipmap,
			DepthBlueNoise,
			DepthFloyd,
			DepthFloydExpo,
			Delta
		};
	}
	typedef SeedModes::Type SeedMode;

	namespace ColorSpaces
	{
		enum Type {
			RGB, HSV, LAB, HN
		};
	}
	typedef ColorSpaces::Type ColorSpace;

	struct Parameters
	{
		Parameters();

		/** camera parameters */
		Camera camera;

		ColorSpace color_space;

		float weight_color;
		float weight_spatial;
		float weight_normal;
		float weight_depth;
		float weight_image;

		/** Number of iterations for superpixel k-means clustering */
		unsigned int iterations;

		/** Superpixel cluster search radius factor */
		float coverage;

		/** Desired radius of a surface element */
		float base_radius;

		/** Desired number of superpixels */
		unsigned int count;

		/** Method used to compute seed points */
		SeedMode seed_mode;

		bool gradient_adaptive_density;
		bool use_density_depth;

		/** Ignores pixels which are too far away or where the depth gradient is too big */
		bool ignore_pixels_with_bad_visibility;

		bool is_conquer_enclaves;

		float segment_threshold;

		bool is_repair_depth;
		bool is_smooth_depth;
		bool is_improve_seeds;

		/** Pixel scala at depth
		 * Radius [px] of a surface element of size base radius [m] and
		 * at given depth [kinect] on the image sensor
		 */
		float computePixelScala(uint16_t depth) const {
			return (depth == 0) ? 0.0f : (camera.focal / camera.convertKinectToMeter(depth) * base_radius);
		}

	};

}

#endif
