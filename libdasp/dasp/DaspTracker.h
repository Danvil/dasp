/*
 * DaspTracker.h
 *
 *  Created on: Feb 14, 2012
 *      Author: david
 */

#ifndef DASPTRACKER_H_
#define DASPTRACKER_H_
//----------------------------------------------------------------------------//
#include "Superpixels.hpp"
#include "SuperpixelHistogram.hpp"
#include "Plots.hpp"
#include <Slimage/Slimage.hpp>
#include <Danvil/Statistics/GMM.h>
#include <boost/thread.hpp>
#include <vector>
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

class DaspTracker
{
public:
	DaspTracker();

	virtual ~DaspTracker();

	/**
	 * Memory contents for raw_kinect_depth and raw_kinect_color must not be
	 * altered until the function step has returned.
	 */
	void step(const slimage::Image1ui16& raw_kinect_depth, const slimage::Image3ub& raw_kinect_color);

	std::map<std::string, slimage::ImagePtr> getImages() const;

	slimage::Image1ub getResultImage() const;

	void Render() const;

	void RenderClusterMap() const;

private:
	void performSegmentationStep();

	void trainInitialColorModel();

	void performTrackingStep();

public:
	bool training_;

	bool enable_smooth_depth_;

	bool show_points_;
	bool show_clusters_;
	bool show_cluster_borders_;
	plots::ColorMode point_color_mode_;
	plots::ColorMode cluster_color_mode_;
	plots::ClusterMode cluster_mode_;
	bool show_graph_;
	bool plot_density_;
	bool plot_segments_;

    boost::shared_ptr<dasp::Parameters> dasp_params;

    float color_model_sigma_scale_;

    unsigned int thread_pool_index_;

private:
	slimage::Image1ui16 kinect_depth;
	slimage::Image3ub kinect_color_rgb;

	slimage::Image1ub result_;

	std::vector<dasp::Seed> seeds;

	Clustering clustering_;

	bool has_hand_gmm_model_;
	boost::shared_ptr<SuperpixelHistogramModel> model_;
	Danvil::GMM::GaussianMixtureModel<2,3,float> hand_gmm_model_;

    plots::ClusterSelection selection_;

	std::map<std::string, slimage::ImagePtr> images_;

	mutable boost::mutex images_mutex_;

	mutable boost::mutex render_mutex_;

};

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
#endif
