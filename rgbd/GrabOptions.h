/*
 * GrabOptions.h
 *
 *  Created on: Aug 22, 2011
 *      Author: david
 */

#ifndef DASP_KINECT_GRABOPTIONS_H_
#define DASP_KINECT_GRABOPTIONS_H_
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

struct GrabOptions
{
	GrabOptions() {
		use_frame_range_ = false;
		use_frame_skip_ = false;
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

	bool use_frame_range_;
	unsigned int frame_first_;
	unsigned int frame_last_;

	bool use_frame_skip_;
	unsigned int frame_skip_;
};

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
