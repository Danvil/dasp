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
//	if(!image_.IsCapabilitySupported(XN_CAPABILITY_ALTERNATIVE_VIEW_POINT)) {
//		cerr << "Can not change image node viewport to depth image!" << endl;
//	}
//	else {
//		image_.GetAlternativeViewPointCap().SetViewPoint(depth_);
//	}
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

	// convert color image to Danvil image
	image_.GetMetaData(imageMD_);
	const XnUInt8* pImage = imageMD_.Data();
	last_color_.resize(imageMD_.FullXRes(), imageMD_.FullYRes());
	memcpy(last_color_.begin(), pImage, last_color_.getByteCount());

	// convert depth image to Danvil image
	depth_.GetMetaData(depthMD_);
	const XnDepthPixel* pDepth = depthMD_.Data();
	last_depth_.resize(depthMD_.FullXRes(), depthMD_.FullYRes());
	memcpy(last_depth_.begin(), pDepth, last_depth_.getByteCount());

//	if(!save_color_fn_.empty()) {
//		Danvil::Images::ImageIO::Save(last_color_, save_color_fn_);
//		save_color_fn_ = "";
//	}
//
//	if(!save_depth_fn_.empty()) {
//		// convert depth to ub1 image for saving as 8bit image
//		Danvil::Images::Image1ubPtr img_save = Danvil::Images::ImageFactory::FactorSimpleImage1ub(depthMD_.FullXRes(), depthMD_.FullYRes());
//		Danvil::memops::convert_uint16_to_uint8<4>(pDepth, img_save->elementCount(), img_save->begin());
//		Danvil::Images::ImageIO::Save(img_save, save_depth_fn_);
//		save_depth_fn_ = "";
//	}

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
