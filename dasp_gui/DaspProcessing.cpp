/*
 * DaspProcessing.cpp
 *
 *  Created on: Feb 14, 2012
 *      Author: david
 */

#include "DaspProcessing.h"
#include <dasp/Neighbourhood.hpp>
#include <dasp/Segmentation.hpp>
#include <dasp/Plots.hpp>
#include <Slimage/Paint.hpp>
#include <Slimage/Convert.hpp>
#define DANVIL_ENABLE_BENCHMARK
#include <Danvil/Tools/Benchmark.h>
#include <Danvil/Color.h>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <stdexcept>
#include <fstream>
using namespace std;
using namespace dasp;
//----------------------------------------------------------------------------//

DaspProcessing::DaspProcessing()
{
	dasp_params.reset(new dasp::Parameters());
	dasp_params->camera = Camera{318.39f, 271.99f, 528.01f, 0.001f};
	dasp_params->seed_mode = dasp::SeedModes::DepthMipmap;
	dasp_params->base_radius = 0.02f;
	dasp_params->gradient_adaptive_density = true;

	color_model_sigma_scale_ = 1.0f;
	thread_pool_index_ = 100;

	show_points_ = false;
	show_clusters_ = true;
	show_cluster_borders_ = true;
	cluster_color_mode_ = plots::Color;
	point_color_mode_ = plots::Color;
	cluster_mode_ = plots::ClusterPoints;
	show_graph_ = false;
	plot_density_ = false;
	plot_segments_ = false;

}

DaspProcessing::~DaspProcessing()
{
}

void DaspProcessing::step(const slimage::Image1ui16& raw_kinect_depth, const slimage::Image3ub& raw_kinect_color)
{
	DANVIL_BENCHMARK_START(step)

	kinect_depth = raw_kinect_depth;
	kinect_color_rgb = raw_kinect_color;

	DANVIL_BENCHMARK_STOP(step)

	performSegmentationStep();

//	DANVIL_BENCHMARK_PRINTALL_COUT
}

slimage::Image3ub ColorizeIntensity(const slimage::Image1f& I, float min, float max, unsigned int pool_id, Danvil::Palette pal=Danvil::Palettes::Blue_Red_Yellow_White)
{
	slimage::Image3ub col(I.width(), I.height());
	if(I) {
		Danvil::ContinuousIntervalColorMapping<unsigned char, float> cm
				= Danvil::ContinuousIntervalColorMapping<unsigned char, float>::Factor(pal);
		cm.useCustomBorderColors(Danvil::Colorub::Black, Danvil::Colorub::White);
		cm.setRange(min, max);
		slimage::ParallelProcess(I, col, [&cm](const slimage::It1f& src, const slimage::It3ub& dst) {
			cm(*src).writeRgb(dst.pointer());
		}, slimage::ThreadingOptions::UsePool(pool_id));
	}
	return col;
}

typedef boost::accumulators::accumulator_set<
	unsigned int,
	boost::accumulators::stats<boost::accumulators::tag::variance>
> CoverageAccType;
CoverageAccType coverage;

void DaspProcessing::performSegmentationStep()
{
	// superpixel parameters
	clustering_.opt = *dasp_params;

	{	boost::interprocess::scoped_lock<boost::mutex> lock(render_mutex_);
		ComputeSuperpixelsIncremental(clustering_, kinect_color_rgb, kinect_depth);
		if(show_clusters_ && (cluster_color_mode_ == plots::CoverageError)) {
			clustering_.ComputeExt();
		}
	}

	DANVIL_BENCHMARK_START(mog)

	slimage::Image1f probability;
	slimage::Image1f probability_2;
	slimage::Image3ub plot_labels;

	result_.resize(kinect_color_rgb.width(), kinect_color_rgb.height());
	result_.fill({0});

	DANVIL_BENCHMARK_STOP(mog)

	DANVIL_BENCHMARK_START(segmentation)
	ClusterLabeling dasp_segment_labeling;
	EdgeWeightGraph dasp_segment_graph;
	if(plot_segments_) {
		// create segmentation graph
		//segments = MinCutSegmentation(clustering_);
		BorderPixelGraph Gnb = CreateNeighborhoodGraph(clustering_, NeighborGraphSettings::SpatialCut());
		EdgeWeightGraph Gnb_local_weights = ComputeEdgeWeights(clustering_, Gnb,
				ClassicSpectralAffinity<true>(clustering_.clusterCount(), clustering_.opt.base_radius));
		dasp_segment_graph = SpectralSegmentation(Gnb_local_weights, boost::get(boost::edge_weight, Gnb_local_weights));
		dasp_segment_labeling = ComputeSegmentLabels(dasp_segment_graph, clustering_.opt.segment_threshold);
	}
	DANVIL_BENCHMARK_STOP(segmentation)

	DANVIL_BENCHMARK_START(plotting)

	{
		slimage::Image3ub vis_img;
		if(show_points_) {
			vis_img = plots::PlotPoints(clustering_, point_color_mode_);
		}
		else {
			vis_img.resize(clustering_.width(), clustering_.height());
			vis_img.fill({{0,0,0}});
		}
		if(show_clusters_) {
			plots::PlotClusters(vis_img, clustering_, cluster_mode_, cluster_color_mode_);
		}
		if(show_cluster_borders_) {
			slimage::Pixel3ub border_color;
			if(cluster_color_mode_ == plots::ColorMode::UniBlack || cluster_color_mode_ == plots::ColorMode::Gradient) {
				border_color = slimage::Pixel3ub{{255,255,255}};
			}
			else if(cluster_color_mode_ == plots::ColorMode::UniWhite || cluster_color_mode_ == plots::ColorMode::Depth) {
				border_color = slimage::Pixel3ub{{0,0,0}};
			}
			else {
				border_color = slimage::Pixel3ub{{255,0,0}};
			}
			plots::PlotEdges(vis_img, clustering_.ComputeLabels(), border_color, 2);
		}

		if(show_graph_ || plot_segments_) {
			if(plot_segments_) {
				// plot segmentation graph
				//std::vector<slimage::Pixel3ub> colors = ComputeSegmentColors(clustering_, labeling);
				std::vector<slimage::Pixel3ub> colors = plots::CreateRandomColors(dasp_segment_labeling.num_labels);
				vis_img = CreateLabelImage(clustering_, dasp_segment_labeling, colors);

				if(show_graph_) {
					plots::PlotWeightedGraphLines(vis_img, clustering_, dasp_segment_graph, [](float weight) {
						return plots::IntensityColor(weight, 0, 1);
					});
				}
			}
			else {
				if(show_graph_) {
					// create neighbourhood graph
					BorderPixelGraph Gnb = CreateNeighborhoodGraph(clustering_);
					// plot neighbourhood graph
					plots::PlotGraphLines(vis_img, clustering_, Gnb);
//					plots::PlotGraphLines(vis_img, clustering_, Gnb, [&Gnb](BorderPixelGraph::edge_descriptor eid) {
//						unsigned int n = boost::get(borderpixels_t(), Gnb, eid).size();
//						return plots::IntensityColor(static_cast<float>(n), 0.0f, 20.0f);
//					});
				}
			}
		}

		// plot label probability
		slimage::Image3ub probability_color;
		slimage::Image3ub probability_2_color;
		if(probability) {
			probability_color = ColorizeIntensity(probability, 0.0f, 1.0f, thread_pool_index_);
		}
		if(probability_2) {
			probability_2_color = ColorizeIntensity(probability_2, 0.0f, 1.0f, thread_pool_index_);
		}

		// visualize density, seed density and density error
		slimage::Image3ub vis_density;
		slimage::Image1ub vis_seed_density;
		slimage::Image3ub vis_density_delta;
		if(plot_density_) {
			slimage::Image1f density = ComputeDepthDensity(clustering_.points, clustering_.opt);
			vis_density = slimage::Image3ub(density.dimensions());
			for(unsigned int i=0; i<vis_density.size(); i++) {
				vis_density[i] = plots::IntensityColor(density[i], 0.0f, 0.1f);
			}
			slimage::Image1f seed_density = ComputeDepthDensityFromSeeds(clustering_.seeds_previous, density);
			slimage::conversion::Convert(seed_density, vis_seed_density); // FIXME * 20.0f
			vis_density_delta.resize(density.width(), density.height());
			for(unsigned int i=0; i<density.size(); i++) {
				vis_density_delta[i] = plots::PlusMinusColor(density[i] - seed_density[i], 0.025f);
			}
		}

		// set images for gui
		{
			boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);

			images_["2D"] = slimage::Ptr(vis_img);
			images_["rhoE"] = slimage::Ptr(vis_density);
			images_["rhoA"] = slimage::Ptr(vis_seed_density);
			images_["rhoD"] = slimage::Ptr(vis_density_delta);
//			images_["seeds"] = slimage::Ptr(seeds_img);
			images_["prob"] = slimage::Ptr(probability_color);
			images_["prob2"] = slimage::Ptr(probability_2_color);
			images_["labels"] = slimage::Ptr(plot_labels);

			for(auto p : sDebugImages) {
				images_[p.first] = p.second;
			}
		}

	}
	DANVIL_BENCHMARK_STOP(plotting)
}

std::map<std::string, slimage::ImagePtr> DaspProcessing::getImages() const
{
	boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);
	return images_;
}

slimage::Image1ub DaspProcessing::getResultImage() const
{
	boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);
	return result_;
}

void DaspProcessing::Render() const
{
	{	boost::interprocess::scoped_lock<boost::mutex> lock(render_mutex_);
		plots::RenderClusters(clustering_, cluster_color_mode_, selection_);
	}
}

void DaspProcessing::RenderClusterMap() const
{
	{	boost::interprocess::scoped_lock<boost::mutex> lock(render_mutex_);
		plots::RenderClusterMap(clustering_, cluster_color_mode_, selection_);
	}
}

//----------------------------------------------------------------------------//
