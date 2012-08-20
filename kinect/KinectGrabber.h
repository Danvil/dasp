/*
 * KinectGrabber.h
 *
 *  Created on: Jun 24, 2011
 *      Author: david
 */

#ifndef ROMOE_KINECT_KINECTGRABBER_H_
#define ROMOE_KINECT_KINECTGRABBER_H_
//----------------------------------------------------------------------------//
#include "GrabOptions.h"
#include <Slimage/Slimage.hpp>
#include <XnCppWrapper.h>
#include <boost/signals.hpp>
#include <string>
namespace xn {
	class Context;
	class DepthGenerator;
	class ImageGenerator;
	class DepthMetaData;
	class ImageMetaData;
}

slimage::Image1ub ColorizeDepth(const slimage::Image1ui16& depth);


//----------------------------------------------------------------------------//
namespace Romeo {
namespace Kinect {
//----------------------------------------------------------------------------//

class KinectGrabber
{
public:
	KinectGrabber();
	virtual ~KinectGrabber();

	void OpenConfig(const std::string& fn_config);

	void OpenFile(const std::string& fn_oni);

	boost::signal<void(slimage::Image3ub)> on_color_;

	boost::signal<void(slimage::Image1ui16)> on_depth_;

	boost::signal<void(slimage::Image1ui16, slimage::Image3ub)> on_depth_and_color_;

	GrabOptions& options() {
		return options_;
	}

	void SetOptions(const GrabOptions& options) {
		options_ = options;
	}

	slimage::Image3ub GetLastColor() const {
		return last_color_;
	}

	slimage::Image1ui16 GetLastDepth() const {
		return last_depth_;
	}

	void GrabAndNotify();

	void Run();

	void Stop() {
		stop_requested_ = true;
	}

private:
	void Init();
	
	bool Grab();

	void Notify();

	bool has_depth_;
	bool has_image_;
	bool can_grab_;

	unsigned int frame_current_;
	GrabOptions options_;

	slimage::Image3ub last_color_;
	slimage::Image1ui16 last_depth_;

	bool stop_requested_;

private:
	xn::Context context_;
	xn::DepthGenerator depth_;
	xn::ImageGenerator image_;
	xn::DepthMetaData depthMD_;
	xn::ImageMetaData imageMD_;
};

//----------------------------------------------------------------------------//
}}
//----------------------------------------------------------------------------//
#endif
