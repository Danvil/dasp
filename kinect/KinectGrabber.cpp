/*
 * KinectGrabber.cpp
 *
 *  Created on: Jun 24, 2011
 *      Author: david
 */

#include "KinectGrabber.h"
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <XnCppWrapper.h>
using namespace std;

slimage::Image1ub ColorizeDepth(const slimage::Image1ui16& depth)
{
	slimage::Image1ub color(depth.dimensions());
	for(unsigned int i=0; i<depth.size(); i++) {
		int v = depth[i];
		int q = ((v - 500)*256) / 2500 + (v % 25);
		unsigned char g = static_cast<unsigned char>(std::min(255, std::max(0, q)));
		color[i] = g;
	}
	return color;
}

//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

struct KinectGrabber::xn_vars
{
	xn::Context context_;
	xn::Player player_;
	xn::DepthGenerator depth_;
	xn::ImageGenerator image_;
	xn::DepthMetaData depthMD_;
	xn::ImageMetaData imageMD_;
};

KinectGrabber::KinectGrabber()
{
	impl_ = std::make_shared<xn_vars>();
	frame_current_ = 0;
	stop_requested_ = false;
}

KinectGrabber::~KinectGrabber()
{
}

void KinectGrabber::OpenConfig(const std::string& fn_config)
{
	xn::EnumerationErrors errors;
	XnStatus rc = impl_->context_.InitFromXmlFile(fn_config.c_str(), &errors);
	if(rc == XN_STATUS_NO_NODE_PRESENT) {
		cerr << "ERROR " << rc << ": " << xnGetStatusString(rc) << endl;
		XnChar strError[1024];
		errors.ToString(strError, 1024);
		cerr << strError << endl;
		throw rc;
	}
	else if(rc != XN_STATUS_OK) {
		cerr << "ERROR " << rc << ": " << xnGetStatusString(rc) << endl;
		throw rc;
	}
	is_oni_ = false;
	Init();
}

void KinectGrabber::OpenFile(const std::string& fn_oni)
{
	XnStatus rc = impl_->context_.Init();
	if(rc != XN_STATUS_OK) {
		cerr << "ERROR " << rc << ": " << xnGetStatusString(rc) << endl;
		throw rc;
	}
	rc = impl_->context_.OpenFileRecording(fn_oni.c_str());
	if(rc != XN_STATUS_OK) {
		cerr << "ERROR " << rc << ": " << xnGetStatusString(rc) << endl;
		throw rc;
	}
	
	// disable looping
	rc = impl_->context_.FindExistingNode(XN_NODE_TYPE_PLAYER, impl_->player_);
	impl_->player_.SetRepeat(false);

	is_oni_ = true;
	Init();
}

void KinectGrabber::Init()
{
	XnStatus rc = impl_->context_.FindExistingNode(XN_NODE_TYPE_DEPTH, impl_->depth_);
	if(rc != XN_STATUS_OK) {
		cout << "Could not get image node. " << rc << ": " << xnGetStatusString(rc) << endl;
		has_depth_ = false;
	}
	else {
		has_depth_ = true;
		impl_->depth_.GetMetaData(impl_->depthMD_);
		cout << "Depth node: " << impl_->depthMD_.FullXRes() << "x" << impl_->depthMD_.FullYRes() << endl;
	}
	rc = impl_->context_.FindExistingNode(XN_NODE_TYPE_IMAGE, impl_->image_);
	if(rc != XN_STATUS_OK) {
		cout << "Could not get image node. " << rc << ": " << xnGetStatusString(rc) << endl;
		has_image_ = false;
	}
	else {
		has_image_ = true;
		impl_->image_.GetMetaData(impl_->imageMD_);
		cout << "Image node: " << impl_->imageMD_.FullXRes() << "x" << impl_->imageMD_.FullYRes() << endl;
	}
	// set color viewport to depth
	if(!impl_->depth_.IsCapabilitySupported(XN_CAPABILITY_ALTERNATIVE_VIEW_POINT)) {
		cerr << "Can not change image node viewport to depth image!" << endl;
	}
	else {
		impl_->depth_.GetAlternativeViewPointCap().SetViewPoint(impl_->image_);
	}
	cout << "Kinect initialized successfully" << endl;
	can_grab_ = true;
}

int KinectGrabber::NumFrames() const
{
	if(is_oni_) {
		XnStatus rc;
		XnUInt32 num1;
		rc = impl_->player_.GetNumFrames(impl_->depth_.GetName(), num1);
		if(rc != XN_STATUS_OK) {
			cerr << "Could not get number of frames for depth node! " << rc << ": " << xnGetStatusString(rc) << endl;
		}
		// XnUInt32 num2;
		// rc = player_.GetNumFrames(image_.GetName(), num2);
		// if(rc != XN_STATUS_OK) {
		// 	cerr << "Could not get number of frames for image node! " << rc << ": " << xnGetStatusString(rc) << endl;
		// }
		// if(num1 != num2) {
		// 	cerr << "Number of frames for depth node (" << num1 << ") and image node (" << num2 << ") are not identical!" << endl;
		// }
		return static_cast<int>(num1);
	}
	else {
		return -1;
	}
}

void KinectGrabber::SeekToFrame(int frame)
{
	if(is_oni_) {
		XnStatus rc;
		// the following line is a hack because calling seektoframe(i) grab seektoframe(i) grab
		// may not yield the same result the second time.
		impl_->player_.SeekToFrame(impl_->depth_.GetName(), (frame == 0 ? 1 : frame - 1), XN_PLAYER_SEEK_SET);
		rc = impl_->player_.SeekToFrame(impl_->depth_.GetName(), frame, XN_PLAYER_SEEK_SET);
		if(rc != XN_STATUS_OK) {
			cerr << "Could not seek to frame in depth node! " << rc << ": " << xnGetStatusString(rc) << endl;
		}
		// rc = player_.SeekToFrame(image_.GetName(), frame, XN_PLAYER_SEEK_SET);
		// if(rc != XN_STATUS_OK) {
		// 	cerr << "Could not seek to frame in image node! " << rc << ": " << xnGetStatusString(rc) << endl;
		// }
	}
	else {
		cerr << "Can not seek in live mode!" << endl;
	}
}

int KinectGrabber::TellFrame()
{
	if(is_oni_) {
		XnStatus rc;
		XnUInt32 num1;
		rc = impl_->player_.TellFrame(impl_->depth_.GetName(), num1);
		if(rc != XN_STATUS_OK) {
			cerr << "Could not tell frame for depth node! " << rc << ": " << xnGetStatusString(rc) << endl;
		}
		// XnUInt32 num2;
		// rc = player_.TellFrame(image_.GetName(), num2);
		// if(rc != XN_STATUS_OK) {
		// 	cerr << "Could not tell frame for image node! " << rc << ": " << xnGetStatusString(rc) << endl;
		// }
		// if(num1 != num2) {
		// 	cerr << "Current frame for depth node (" << num1 << ") and image node (" << num2 << ") are not identical!" << endl;
		// }
		return static_cast<int>(num1);
	}
	else {
		return -1;
	}
}

bool KinectGrabber::GrabAndNotify()
{
	if(can_grab_) {
		can_grab_ = Grab();
	}
	if(can_grab_) {
		Notify();
	}
	return can_grab_;
}

void KinectGrabber::Run()
{
	while(!stop_requested_ && can_grab_) {
		GrabAndNotify();
	}
}

bool KinectGrabber::Grab()
{
	// read a new frame
	XnStatus rc = impl_->context_.WaitAndUpdateAll();
	if(rc != XN_STATUS_OK) {
		cerr << "ERROR " << rc << ": " << xnGetStatusString(rc) << endl;
		return false;
	}

	// only process frames in frame range
	if(options_.use_frame_range_) {
		if(frame_current_ < options_.frame_first_) {
			frame_current_++;
			return true;
		}
		if(options_.frame_last_ < frame_current_) {
			return false;
		}
	}

	// skip frames if desired by user
	if(options_.use_frame_skip_ && frame_current_ % options_.frame_skip_ != 0) {
		frame_current_++;
		return true;
	}

	// store color image
	impl_->image_.GetMetaData(impl_->imageMD_);
	const XnUInt8* pImage = impl_->imageMD_.Data();
	last_color_.resize(impl_->imageMD_.FullXRes(), impl_->imageMD_.FullYRes());
	memcpy(last_color_.begin().pointer(), pImage, last_color_.size() * 3 * sizeof(unsigned char));

	// store depth image
	impl_->depth_.GetMetaData(impl_->depthMD_);
	const XnDepthPixel* pDepth = impl_->depthMD_.Data();
	last_depth_.resize(impl_->depthMD_.FullXRes(), impl_->depthMD_.FullYRes());
	memcpy(last_depth_.begin().pointer(), pDepth, last_depth_.size() * sizeof(uint16_t));

	frame_current_++;
	return true;
}

void KinectGrabber::Notify()
{
	if(!on_color_.empty()) {
		on_color_(last_color_);
	}
	if(!on_depth_.empty()) {
		on_depth_(last_depth_);
	}
	if(!on_depth_and_color_.empty()) {
		on_depth_and_color_(last_depth_, last_color_);
	}
}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
