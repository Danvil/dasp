/*
 * GrabOptions.h
 *
 *  Created on: Aug 22, 2011
 *      Author: david
 */

#ifndef ROMEO_KINECT_GRABOPTIONS_H_
#define ROMEO_KINECT_GRABOPTIONS_H_
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
namespace Romeo {
namespace Kinect {
//----------------------------------------------------------------------------//

struct GrabOptions
{
	GrabOptions() {
		use_frame_range_ = false;
		use_frame_skip_ = false;
		use_z_range_ = false;
	}

	void EnableFrameRange(unsigned int first, unsigned int last) {
		use_frame_range_ = true;
		frame_first_ = first;
		frame_last_ = last;
	}

	void EnableFrameSkip(unsigned int skip) {
		use_frame_skip_ = true;
		frame_skip_ = skip;
	}

	void EnableDepthRange(float z_min, float z_max) {
		use_z_range_ = true;
		z_min_ = z_min;
		z_max_ = z_max;
	}

	bool use_frame_range_;
	unsigned int frame_first_;
	unsigned int frame_last_;

	bool use_frame_skip_;
	unsigned int frame_skip_;

	bool use_z_range_;
	float z_min_;
	float z_max_;
};

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
#endif
