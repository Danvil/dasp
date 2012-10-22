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
namespace Romeo {
namespace Kinect {
//----------------------------------------------------------------------------//

KinectGrabber::KinectGrabber()
{
	frame_current_ = 0;
	stop_requested_ = false;
}

KinectGrabber::~KinectGrabber()
{
}

void KinectGrabber::OpenConfig(const std::string& fn_config)
{
	xn::EnumerationErrors errors;
	XnStatus rc = context_.InitFromXmlFile(fn_config.c_str(), &errors);
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
	XnStatus rc = context_.Init();
	if(rc != XN_STATUS_OK) {
		cerr << "ERROR " << rc << ": " << xnGetStatusString(rc) << endl;
		throw rc;
	}
	rc = context_.OpenFileRecording(fn_oni.c_str());
	if(rc != XN_STATUS_OK) {
		cerr << "ERROR " << rc << ": " << xnGetStatusString(rc) << endl;
		throw rc;
	}
	
	// disable looping
	rc = context_.FindExistingNode(XN_NODE_TYPE_PLAYER, player_);
	player_.SetRepeat(false);

	is_oni_ = true;
	Init();
}

void KinectGrabber::Init()
{
	XnStatus rc = context_.FindExistingNode(XN_NODE_TYPE_DEPTH, depth_);
	if(rc != XN_STATUS_OK) {
		cout << "Could not get image node. " << rc << ": " << xnGetStatusString(rc) << endl;
		has_depth_ = false;
	}
	else {
		has_depth_ = true;
		depth_.GetMetaData(depthMD_);
		cout << "Depth node: " << depthMD_.FullXRes() << "x" << depthMD_.FullYRes() << endl;
	}
	rc = context_.FindExistingNode(XN_NODE_TYPE_IMAGE, image_);
	if(rc != XN_STATUS_OK) {
		cout << "Could not get image node. " << rc << ": " << xnGetStatusString(rc) << endl;
		has_image_ = false;
	}
	else {
		has_image_ = true;
		image_.GetMetaData(imageMD_);
		cout << "Image node: " << imageMD_.FullXRes() << "x" << imageMD_.FullYRes() << endl;
	}
	// set color viewport to depth
	if(!depth_.IsCapabilitySupported(XN_CAPABILITY_ALTERNATIVE_VIEW_POINT)) {
		cerr << "Can not change image node viewport to depth image!" << endl;
	}
	else {
		depth_.GetAlternativeViewPointCap().SetViewPoint(image_);
	}
	cout << "Kinect initialized successfully" << endl;
	can_grab_ = true;
}

int KinectGrabber::NumFrames() const
{
	if(is_oni_) {
		XnStatus rc;
		XnUInt32 num1;
		rc = player_.GetNumFrames(depth_.GetName(), num1);
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
		player_.SeekToFrame(depth_.GetName(), (frame == 0 ? 1 : frame - 1), XN_PLAYER_SEEK_SET);
		rc = player_.SeekToFrame(depth_.GetName(), frame, XN_PLAYER_SEEK_SET);
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
		rc = player_.TellFrame(depth_.GetName(), num1);
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
	XnStatus rc = context_.WaitAndUpdateAll();
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
	image_.GetMetaData(imageMD_);
	const XnUInt8* pImage = imageMD_.Data();
	last_color_.resize(imageMD_.FullXRes(), imageMD_.FullYRes());
	memcpy(last_color_.begin().pointer(), pImage, last_color_.size() * 3 * sizeof(unsigned char));

	// store depth image
	depth_.GetMetaData(depthMD_);
	const XnDepthPixel* pDepth = depthMD_.Data();
	last_depth_.resize(depthMD_.FullXRes(), depthMD_.FullYRes());
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
}}
//----------------------------------------------------------------------------//
