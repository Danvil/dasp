/*
 * DaspTracker.h
 *
 *  Created on: Feb 14, 2012
 *      Author: david
 */

#ifndef DASPTRACKER_H_
#define DASPTRACKER_H_
//----------------------------------------------------------------------------//
#include <SuperPoints/Superpixels.hpp>
#include <SuperPoints/SuperpixelHistogram.hpp>
#include <Slimage/Slimage.hpp>
#include <Danvil/Statistics/GMM.h>
#include <Danvil/Images/Image.h>
#include <boost/thread.hpp>
#include <vector>
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

class DaspTracker
{
public:
	DaspTracker();

	virtual ~DaspTracker();

	void step(Danvil::Images::Image1ui16Ptr raw_kinect_depth, Danvil::Images::Image3ubPtr raw_kinect_color);

	std::map<std::string, slimage::ImagePtr> getImages() const;

private:
	void performSegmentationStep();

	void trainInitialColorModel();

	void performTrackingStep();

public:
	bool training_;

    boost::shared_ptr<dasp::Parameters> dasp_params;

    float color_model_sigma_scale_;

private:
	slimage::Image1ui16 kinect_depth;
	slimage::Image3ub kinect_color_rgb;

	slimage::Image3f kinect_color;

	dasp::ImagePoints points;
	std::vector<dasp::Cluster> clusters;

	bool has_hand_gmm_model_;
	Danvil::GMM::GaussianMixtureModel<5,3,float> hand_gmm_model_;

	std::vector<dasp::SuperpixelHistogram> model_hist_;

	std::map<std::string, slimage::ImagePtr> images_;

	mutable boost::mutex images_mutex_;

};

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
#endif
