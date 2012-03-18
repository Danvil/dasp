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

	void update(double dt, double current_time) {
		GrabAndNotify();
	}

	bool Grab();

	void Run();

	boost::signal<void(slimage::Image3ub)> on_color_;

	boost::signal<void(slimage::Image1ui16)> on_depth_;

	boost::signal<void(slimage::Image1ui16, slimage::Image3ub)> on_depth_and_color_;

	std::string save_color_fn_;
	std::string save_depth_fn_;
	std::string save_points_fn_;

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

	bool stop_requested_;

private:
	void GrabAndNotify();

	void Notify();

private:
	void Init();
	bool has_depth_;
	bool has_image_;
	bool can_grab_;

	unsigned int frame_current_;
	GrabOptions options_;

	slimage::Image3ub last_color_;
	slimage::Image1ui16 last_depth_;

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
