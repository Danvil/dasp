/*
 * KinectGrabberObject.hpp
 *
 *  Created on: Apr 17, 2012
 *      Author: david
 */

#ifndef KINECTGRABBEROBJECT_HPP_
#define KINECTGRABBEROBJECT_HPP_

#include <KinectGrabber.h>
#include <Danvil/SimpleEngine.h>
#include <boost/bind.hpp>

class KinectGrabberObject
: public Danvil::SimpleEngine::IUpdateable
{
public:
	KinectGrabberObject() {
		kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());

		std::cout << "Initializing kinect..." << std::endl;
		kinect_grabber_->OpenConfig("/home/david/Programs/RGBD/OpenNI/Platform/Linux-x86/Redist/Samples/Config/SamplesConfig.xml");

//		kinect_grabber_->on_depth_and_color_.connect(boost::bind(&KinectGrabberObject::notify, this, _1, _2));
	}

	void setNotification(const boost::function<void(slimage::Image1ui16, slimage::Image3ub)>& f) {
		fnc_ = f;
	}

//	void notify(slimage::Image1ui16 depth, slimage::Image3ub color) {
//		fnc_(depth, color);
//	}

	void update(double dt, double time) {
//		kinect_grabber_->GrabAndNotify();
		for(unsigned int i=0; i<10; i++) {
			kinect_grabber_->Grab();
		}
		if(fnc_) {
			fnc_(kinect_grabber_->GetLastDepth(), kinect_grabber_->GetLastColor());
		}
	}

private:
	boost::shared_ptr<Romeo::Kinect::KinectGrabber> kinect_grabber_;

	boost::function<void(slimage::Image1ui16, slimage::Image3ub)> fnc_;
};

#endif
