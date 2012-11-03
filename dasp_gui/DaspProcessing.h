/*
 * DaspProcessing.h
 *
 *  Created on: Feb 14, 2012
 *      Author: david
 */

#ifndef DASPGUI_DASPPROCESSING_H_
#define DASPGUI_DASPPROCESSING_H_
//----------------------------------------------------------------------------//
#include <dasp/Superpixels.hpp>
#include <dasp/Plots.hpp>
#include <dasp/Graph.hpp>
#include <Slimage/Slimage.hpp>
#include <boost/thread.hpp>
#include <vector>
//----------------------------------------------------------------------------//

class DaspProcessing
{
public:
	DaspProcessing();

	virtual ~DaspProcessing();

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

public:
	bool show_points_;
	bool show_clusters_;
	bool show_cluster_borders_;
	dasp::plots::ColorMode point_color_mode_;
	dasp::plots::ColorMode cluster_color_mode_;
	dasp::plots::ClusterMode cluster_mode_;
	bool graph_cut_spatial_;
	bool show_graph_;
	bool show_graph_weights_;
	bool plot_segments_;
	bool plot_density_;

    boost::shared_ptr<dasp::Parameters> dasp_params;

    float color_model_sigma_scale_;

    unsigned int thread_pool_index_;

private:
	slimage::Image1ui16 kinect_depth;
	slimage::Image3ub kinect_color_rgb;

	slimage::Image1ub result_;

	std::vector<dasp::Seed> seeds;

public:
	dasp::Superpixels clustering_;

	dasp::BorderPixelGraph Gnb;
	dasp::EdgeWeightGraph Gnb_weighted;

private:
	dasp:: plots::ClusterSelection selection_;

	std::map<std::string, slimage::ImagePtr> images_;

	mutable boost::mutex images_mutex_;

	mutable boost::mutex render_mutex_;

};

//----------------------------------------------------------------------------//
#endif
