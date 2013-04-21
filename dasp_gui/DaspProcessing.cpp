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
#include <dasp/Metric.hpp>
#include <dasp/impl/Sampling.hpp>
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
	dasp_params->seed_mode = dasp::SeedModes::SimplifiedPDS;
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
	graph_cut_spatial_ = true;
	show_graph_ = false;
	show_graph_weights_ = 2;
	plot_segments_ = false;
	plot_density_ = false;

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

	// slimage::Image1f probability;
	// slimage::Image1f probability_2;
	// slimage::Image3ub plot_labels;

	result_.resize(kinect_color_rgb.width(), kinect_color_rgb.height());
	result_.fill({0});

	DANVIL_BENCHMARK_STOP(mog)

	DANVIL_BENCHMARK_START(graph)
	if(show_graph_ || plot_segments_) {
		// create neighbourhood graph
		Gnb = CreateNeighborhoodGraph(clustering_,
			graph_cut_spatial_ ? NeighborGraphSettings::SpatialCut() : NeighborGraphSettings::NoCut());
		Gnb_weighted = ComputeEdgeWeights(clustering_, Gnb,
				DepthAdaptiveMetric(
					clustering_.opt.weight_spatial, clustering_.opt.weight_color, clustering_.opt.weight_normal,
					clustering_.opt.base_radius));
	}
	DANVIL_BENCHMARK_STOP(graph)

	DANVIL_BENCHMARK_START(segmentation)
	UndirectedWeightedGraph similarity_graph;
	UndirectedWeightedGraph dasp_segment_graph;
	graphseg::GraphLabeling dasp_segment_labeling;
	if(plot_segments_ || show_graph_weights_ >= 3) {
		// create segmentation graph
		//segments = MinCutSegmentation(clustering_);
		similarity_graph = ComputeEdgeWeights(clustering_, Gnb,
				ClassicSpectralAffinity<true>(
					clustering_.clusterCount(), clustering_.opt.base_radius,
					1.0f, 2.0f, 3.0f));
		if(plot_segments_ || show_graph_weights_ >= 4) {
//					clustering_.opt.weight_spatial, clustering_.opt.weight_color, clustering_.opt.weight_normal));
			dasp_segment_graph = SpectralSegmentation(similarity_graph, boost::get(boost::edge_bundle, similarity_graph));
			dasp_segment_labeling = graphseg::ComputeSegmentLabels(dasp_segment_graph, clustering_.opt.segment_threshold);
		}
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

		if(plot_segments_) {
			// plot segments
			//std::vector<slimage::Pixel3ub> colors = ComputeSegmentColors(clustering_, labeling);
			std::vector<slimage::Pixel3ub> colors = plots::CreateRandomColors(dasp_segment_labeling.num_labels);
			vis_img = CreateLabelImage(clustering_, dasp_segment_labeling, colors);
		}

		if(show_graph_) {
			switch(show_graph_weights_) {
				default:
				case 0:
				case 1: // FIXME implement
					plots::PlotGraphLines(vis_img, clustering_, Gnb);
					break;
				case 2: // dasp metric
					plots::PlotWeightedGraphLines(vis_img, clustering_, Gnb_weighted,
						[](float weight) {
							return plots::DistanceColor(weight, 0, 10);
						});
					break;
				case 3: // spectral metric
					plots::PlotWeightedGraphLines(vis_img, clustering_, similarity_graph,
						[](float sim) {
							return plots::IntensityColor(sim, 0, 1);
						});
					break;
				case 4: { // spectral result
					const float T = clustering_.opt.segment_threshold;
					plots::PlotWeightedGraphLines(vis_img, clustering_, dasp_segment_graph,
						[T](float sim) {
							float q = std::max(0.0f, 2.0f*T - sim);
							return plots::IntensityColor(q, 0, 2.0f*T);
						});
					break;
				}
				case 5: { // spectral result
					vis_img.fill({{255,255,255}});
					const float T = clustering_.opt.segment_threshold;
					plots::PlotWeightedGraphLines(vis_img, clustering_, dasp_segment_graph,
						[T](float sim) {
							return plots::IntensityColorBW(2.0f*T - sim, 0, 2.0f*T);
						});
					for(int i=0; i<clustering_.points.size(); i++) {
						if(!clustering_.points[i].is_valid) {
							vis_img[i] = slimage::Pixel3ub{{0,0,0}};
						}
					}
					break;
				}
				case 6: { // ucm
					vis_img.fill({{255,255,255}});
					const float T = clustering_.opt.segment_threshold;
					for(auto eid : as_range(boost::edges(dasp_segment_graph))) {
						int i = boost::source(eid, dasp_segment_graph);
						int j = boost::target(eid, dasp_segment_graph);
						float q = dasp_segment_graph[eid] /(2.0f*T);
						unsigned char g = static_cast<unsigned char>(255.0f*(1.0f - q));
						slimage::Pixel3ub color{{g,g,g}};
						dasp::NeighbourhoodGraph::edge_descriptor e = boost::edge(i, j, Gnb).first;
						for(int k : Gnb[e].border_pixel_ids) { 
							vis_img[k] = color;
						}
					}
					for(int i=0; i<clustering_.points.size(); i++) {
						if(!clustering_.points[i].is_valid) {
							vis_img[i] = slimage::Pixel3ub{{0,0,0}};
						}
					}
					break;
				}
			}
		}

		// // plot label probability
		// slimage::Image3ub probability_color;
		// slimage::Image3ub probability_2_color;
		// if(probability) {
		// 	probability_color = ColorizeIntensity(probability, 0.0f, 1.0f, thread_pool_index_);
		// }
		// if(probability_2) {
		// 	probability_2_color = ColorizeIntensity(probability_2, 0.0f, 1.0f, thread_pool_index_);
		// }

		// visualize density, seed density and density error
		slimage::Image3ub vis_density;
		slimage::Image3ub vis_saliency;
		slimage::Image3ub vis_seed_density;
		slimage::Image3ub vis_density_delta;
		if(plot_density_) {
			Eigen::MatrixXf density = ComputeDepthDensity(clustering_.points, clustering_.opt);
			vis_density = slimage::Image3ub(density.rows(), density.cols());
			for(unsigned int i=0; i<vis_density.size(); i++) {
				float d = density.data()[i];
				vis_density[i] = plots::IntensityColor(d, 0.0f, 0.015f);
			}

			Eigen::MatrixXf saliency = clustering_.saliency;
			vis_saliency = slimage::Image3ub(saliency.rows(), saliency.cols());
			for(unsigned int i=0; i<vis_saliency.size(); i++) {
				float d = saliency.data()[i];
				vis_saliency[i] = plots::PlusMinusColor(d, +1.0f);
			}

			// FIXME plot combined density

			if(clustering_.opt.seed_mode == SeedModes::Delta) {
				Eigen::MatrixXf seed_density = ComputeDepthDensityFromSeeds(clustering_.seeds_previous, density);
				vis_seed_density.resize(density.rows(), density.cols());
				for(unsigned int i=0; i<seed_density.size(); i++) {
					//vis_seed_density[i] = static_cast<unsigned char>(255.0f * 20.0f * seed_density[i]);
					float d = seed_density.data()[i];
					vis_seed_density[i] = plots::IntensityColor(d, 0.0f, 0.015f);
				}
//				slimage::conversion::Convert(seed_density, vis_seed_density);

				vis_density_delta.resize(density.rows(), density.cols());
				for(unsigned int i=0; i<density.size(); i++) {
					float q = density.data()[i] - seed_density.data()[i];
					vis_density_delta[i] = plots::PlusMinusColor(q, 0.010f);
				}
			}
		}

		// set images for gui
		{
			boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);

			if(vis_img) images_["2D"] = slimage::Ptr(vis_img);
			if(vis_density) images_["density"] = slimage::Ptr(vis_density);
			if(vis_saliency) images_["saliency"] = slimage::Ptr(vis_saliency);
			if(vis_seed_density) images_["density (seeds)"] = slimage::Ptr(vis_seed_density);
			if(vis_density_delta) images_["density (delta)"] = slimage::Ptr(vis_density_delta);
//			images_["seeds"] = slimage::Ptr(seeds_img);
			// images_["prob"] = slimage::Ptr(probability_color);
			// images_["prob2"] = slimage::Ptr(probability_2_color);
			// images_["labels"] = slimage::Ptr(plot_labels);

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
