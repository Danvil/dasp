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

void KinectGrabber::GrabAndNotify()
{
	if(can_grab_) {
		can_grab_ = Grab();
		Notify();
	}
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
	XnStatus rc = context_.WaitAnyUpdateAll();
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
