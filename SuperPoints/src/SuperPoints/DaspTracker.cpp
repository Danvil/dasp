/*
 * DaspTracker.cpp
 *
 *  Created on: Feb 14, 2012
 *      Author: david
 */

#include "DaspTracker.h"
#include "Mipmaps.hpp"
#include "BlueNoise.hpp"
#include "PointsAndNormals.hpp"
#include "AutoDepth.hpp"
#include "Plots.hpp"
#include <Slimage/Paint.hpp>
#define DANVIL_ENABLE_BENCHMARK
#include <Danvil/Tools/Benchmark.h>
#include <Danvil/Color.h>
#include <Danvil/Color/LAB.h>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <stdexcept>
#include <fstream>
using namespace std;
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

DaspTracker::DaspTracker()
{
	training_ = false;

	dasp_params.reset(new dasp::Parameters());
//	// opencv calibration
//	528.293477 0.000000 318.394546
//	0.000000 527.719472 271.987945
//	0.000000 0.000000 1.000000
	dasp_params->camera = Camera{318.39f, 271.99f, 528.01f, 0.001f};
	dasp_params->seed_mode = dasp::SeedModes::DepthMipmap;
	dasp_params->base_radius = 0.02f;
	dasp_params->gradient_adaptive_density = true;

	color_model_sigma_scale_ = 1.0f;
	thread_pool_index_ = 100;
	has_hand_gmm_model_ = false;

	show_points_ = false;
	show_clusters_ = true;
	show_cluster_borders_ = true;
	cluster_color_mode_ = plots::Color;
	point_color_mode_ = plots::Color;
	cluster_mode_ = plots::ClusterPoints;
	show_graph_ = false;
	plot_density_ = false;

}

DaspTracker::~DaspTracker()
{
}

//void WdgtKinectSuperPoints::ComputeBlueNoiseImpl()
//{
//	slimage::Image1ub img_raw = slimage::Load1ub("/home/david/Documents/DataSets/2012-02-06 Blue Noise/flower_256.png");
//	slimage::Image1f density(img_raw.width(), img_raw.height());
//	for(unsigned int i=0; i<density.size(); i++) {
//		density[i] = 0.25f*(1.0f - float(img_raw[i]) / 255.0f);
//	}
//
//	while(running_) {
//		std::vector<dasp::BlueNoise::Point> points = dasp::BlueNoise::Compute(density);
//		slimage::Image1ub img_pnts(density.width(), density.height());
//		img_pnts.fill(255);
//		dasp::BlueNoise::PlotPoints(points, img_pnts);
//
//		// set images for gui
//		images_mutex_.lock();
//		images_["raw"] = slimage::Ptr(img_raw);
//		images_["blue"] = slimage::Ptr(img_pnts);
//		images_mutex_.unlock();
//	}
//}

void DaspTracker::step(const slimage::Image1ui16& raw_kinect_depth, const slimage::Image3ub& raw_kinect_color)
{
	DANVIL_BENCHMARK_START(step)

	kinect_depth = raw_kinect_depth;

	kinect_color_rgb = raw_kinect_color;

	// convert rgb to lab
	kinect_color.resize(kinect_color_rgb.width(), kinect_color_rgb.height());
	slimage::ParallelProcess(kinect_color_rgb, kinect_color, [](const unsigned char* rgb, float* lab) {
		float r = float(rgb[0]) / 255.0f;
		float g = float(rgb[1]) / 255.0f;
		float b = float(rgb[2]) / 255.0f;

//		Danvil::color_rgb_to_lab(r, g, b, r, g, b);
//		r /= 1000.0f;
//		g /= 100.0f;
//		b /= 100.0f;

//		float a = r + g + b;
//		if(a > 0.05f) {
//			r /= a;
//			g /= a;
//			b /= a;
//		}
//		else {
//			r = 0;
//			g = 0;
//			b = 0;
//		}

//		float a = r + g + b;
//		if(a > 0.05f) {
//			r /= a;
//			g /= a;
//			b = a * 0.1;
//		}
//		else {
//			r = 0;
//			g = 0;
//			b = 0;
//		}

		lab[0] = r;
		lab[1] = g;
		lab[2] = b;
	}, slimage::ThreadingOptions::Single());

	DANVIL_BENCHMARK_STOP(step)

	performSegmentationStep();

	if(training_) {
		trainInitialColorModel();
		training_ = false;
	}

	performTrackingStep();

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
		slimage::ParallelProcess(I, col, [&cm](const float* src, unsigned char* dst) {
			cm(*src).writeRgb(dst);
		}, slimage::ThreadingOptions::UsePool(pool_id));
	}
	return col;
}

typedef boost::accumulators::accumulator_set<
	unsigned int,
	boost::accumulators::stats<boost::accumulators::tag::variance>
> CoverageAccType;
CoverageAccType coverage;

void DaspTracker::performSegmentationStep()
{
	// superpixel parameters
	clustering_.opt = *dasp_params;

	// compute normals only if necessary
	slimage::Image3f kinect_normals;
//	if(super_params_ext.weight_normal > 0.0f) {
//		DANVIL_BENCHMARK_START(normals)
//		slimage::Image3f kinect_points = dasp::PointsAndNormals::ComputePoints(kinect_depth, slimage::ThreadingOptions::UsePool(thread_pool_index_));
//		slimage::Image3f kinect_normals = dasp::PointsAndNormals::ComputeNormals(kinect_depth, kinect_points, slimage::ThreadingOptions::UsePool(thread_pool_index_));
//	//	dasp::PointsAndNormals::ComputeNormalsFast(kinect_depth, kinect_points, kinect_normals);
//		DANVIL_BENCHMARK_STOP(normals)
//	}

	// prepare super pixel points
	DANVIL_BENCHMARK_START(points)
	ImagePoints old_points = clustering_.points;
	clustering_.CreatePoints(kinect_color, kinect_depth, kinect_normals);
	DANVIL_BENCHMARK_STOP(points)

	// compute super pixel seeds
	DANVIL_BENCHMARK_START(seeds)
	std::vector<Seed> old_seeds = clustering_.getClusterCentersAsSeeds();
	seeds = clustering_.FindSeeds(old_seeds, old_points);
//	std::cout << "Seeds: " << seeds.size() << std::endl;
	DANVIL_BENCHMARK_STOP(seeds)

	// compute super pixel point edges and improve seeds with it
//	slimage::Image1f edges;
//	DANVIL_BENCHMARK_START(improve)
//	dasp::ComputeEdges(points, edges, super_params_ext);
//	if(!edges.isNull()) {
//		dasp::ImproveSeeds(seeds, points, edges, super_params_ext);
//	}
//	DANVIL_BENCHMARK_STOP(improve)

	// compute clusters
	DANVIL_BENCHMARK_START(clusters)
	{	boost::interprocess::scoped_lock<boost::mutex> lock(render_mutex_);
		clustering_.ComputeSuperpixels(seeds);
		if(show_clusters_ && (cluster_color_mode_ == plots::Eccentricity || cluster_color_mode_ == plots::Circularity || cluster_color_mode_ == plots::Thickness)) {
			clustering_.ComputeClusterInfo();
		}
	}
	DANVIL_BENCHMARK_STOP(clusters)

	DANVIL_BENCHMARK_START(segmentation)

	slimage::Image1f probability;
	slimage::Image1f probability_2;
	slimage::Image3ub plot_labels;

	if(has_hand_gmm_model_) {

		auto scaled_model = hand_gmm_model_;
		scaled_model.ScaleDeviation(color_model_sigma_scale_); // FIXME get parameter from Gui

		// classify color with gmm
		std::vector<float> cluster_probability = clustering_.ForClusterCenters(
				[&scaled_model](const dasp::Point& center) {
					Danvil::ctLinAlg::Vec3f v(center.color[0], center.color[1], center.color[2]);
					return scaled_model(v);
//					// detect red
//					return v.x / (v.x + v.y + v.z);
				});

		unsigned int cnt_over_active = 0;
		for(float& v : cluster_probability) {
			if(v > 0.50) {
				cnt_over_active ++;
				v = 1.0f;
			}
			else {
				v = 0.0f;
			}
		}
		coverage(cnt_over_active);
		std::cout << "Active cluster: " << cnt_over_active << std::endl;

		{	boost::interprocess::scoped_lock<boost::mutex> lock(render_mutex_);
			selection_ = plots::ClusterSelection::Empty(clustering_.clusterCount());
			for(unsigned int i=0; i<cluster_probability.size(); i++) {
				selection_[i] = (cluster_probability[i] > 0.80f);
			}
		}

		// paint cluster probabilities
		probability.resize(clustering_.width(), clustering_.height());
		probability.fill(0.0f);
		clustering_.ForPixelClusters([&probability,&cluster_probability](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
			probability(p.spatial_x(), p.spatial_y()) = cluster_probability[cid];
		});

//		// histogram color model
//		SuperpixelGraph G = clustering_.CreateNeighborhoodGraph();
//		std::vector<float> cluster_probability_2 = model_->evaluate(G);
//
//		probability_2.resize(clustering_.width(), clustering_.height());
//		probability_2.fill(-1.0f);
//		clustering_.ForPixelClusters([&probability_2,&cluster_probability_2](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
//			float p_base = cluster_probability_2[cid];
//			probability_2(p.spatial_x(), p.spatial_y()) = p_base;
//		});
//
//		// plot model gaussian cluster id for each superpixel
//		std::vector<unsigned int> labels = model_->label(G);
//		plot_labels.resize(kinect_color.width(), kinect_color.height());
//		plot_labels.fill(0);
//		std::vector<Eigen::Vector3f> model_colors = model_->getClusterColors();
//		clustering_.ForPixelClusters([&plot_labels,&labels,&model_colors](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
//			unsigned int l = labels[cid];
//			Eigen::Vector3f color = model_colors[l];
//			plot_labels(p.spatial_x(), p.spatial_y()) = slimage::Pixel3ub{
//				{(unsigned char)(255.0f*color[0]), (unsigned char)(255.0f*color[1]), (unsigned char)(255.0f*color[2])}
//			};
//		});

	}

	result_.resize(kinect_color.width(), kinect_color.height());
	if(has_hand_gmm_model_) {
		result_ = slimage::Convert_f_2_ub(probability);
	}
	else {
		result_.fill(0);
	}

	DANVIL_BENCHMARK_STOP(segmentation)

	DANVIL_BENCHMARK_START(plotting)

	{
//		// plot super pixel color (with edges)
//		slimage::Image3ub vis_super_color = plots::PlotClusters(clustering_, plots::ClusterPoints, plots::Color);
//		std::vector<int> superpixel_labels = clustering_.ComputePixelLabels();
//		plots::PlotEdges(vis_super_color, superpixel_labels, slimage::Pixel3ub{{0, 0, 0}}, 2);
//
//		boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);
//		images_["s_rgb"] = slimage::Ptr(vis_super_color);

//		std::vector<ClusterInfo> c_info = clustering_.ComputeClusterInfo();
//		ClusterGroupInfo cg_info = clustering_.ComputeClusterGroupInfo(c_info);
//		std::cout << "hist_eccentricity: " << cg_info.hist_eccentricity << std::endl;
//		std::cout << "hist_radius: " << cg_info.hist_radius << std::endl;
//		std::cout << "hist_thickness: " << cg_info.hist_thickness << std::endl;
//
//		slimage::Image3ub vis_cluster_circularity(clustering_.width(), clustering_.height());
//		vis_cluster_circularity.fill(0);
//		clustering_.ForPixelClusters([&vis_cluster_circularity,&c_info](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
//			vis_cluster_circularity(p.spatial_x(), p.spatial_y()) = plots::IntensityColor(c_info[cid].coverage, 0.0f, 1.0f);
//		});
//
//		slimage::Image3ub vis_cluster_thickness(clustering_.width(), clustering_.height());
//		vis_cluster_thickness.fill(0);
//		clustering_.ForPixelClusters([&vis_cluster_thickness,&c_info](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
//			vis_cluster_thickness(p.spatial_x(), p.spatial_y()) = plots::IntensityColor(c_info[cid].t, 0.0f, 0.01f);
//		});

		slimage::Image3ub vis_img;
		if(show_points_) {
			vis_img = plots::PlotPoints(clustering_, point_color_mode_);
		}
		else {
			vis_img.resize(clustering_.width(), clustering_.height());
			vis_img.fill(0);
		}
		if(show_clusters_) {
			plots::PlotClusters(vis_img, clustering_, cluster_mode_, cluster_color_mode_);
		}
		if(show_cluster_borders_) {
			plots::PlotEdges(vis_img, clustering_.ComputePixelLabels(), slimage::Pixel3ub{{0,0,0}}, 2);
		}

		// create and plot neighbourhood graph
		if(show_graph_) {
			SuperpixelGraph G = clustering_.CreateNeighborhoodGraph();
			for(unsigned int i=0; i<G.nodes_.size(); i++) {
				for(unsigned int k : G.node_connections_[i]) {
					slimage::PaintLine(vis_img, G.nodes_[i].x, G.nodes_[i].y, G.nodes_[k].x, G.nodes_[k].y, slimage::Pixel3ub{{255,255,255}});
				}
			}
		}

//		// plot point depth
//		slimage::Image3ub vis_point_depth = plots::PlotPoints(clustering_, plots::Depth);
//
//		// plot point normals
//		slimage::Image3ub vis_point_normal = plots::PlotPoints(clustering_, plots::Gradient);
//
////		// plot density and seeds
////		slimage::Image1f density = dasp::ComputeDepthDensity(points, super_params_ext);
////		slimage::Image3ub seeds_img = ColorizeIntensity(density, 0.0f, 0.02f, thread_pool_index_, Danvil::Palettes::Blue_Red_Yellow);
////		dasp::PlotSeeds(seeds, seeds_img, slimage::Pixel3ub{{255,255,255}}, 3);
//
//		// plot superpixel depth
//		slimage::Image3ub vis_super_depth = plots::PlotClusters(clustering_, plots::ClusterPoints, plots::Depth);
//
//		// plot superpixel normals
//		slimage::Image3ub vis_super_normal = plots::PlotClusters(clustering_, plots::ClusterPoints, plots::Gradient);
//
//		// plot superpixel crosses
//		slimage::Image3ub vis_super_cross = plots::PlotClusters(clustering_, plots::ClusterEllipses, plots::Color);

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
		slimage::Image1ub vis_density;
		slimage::Image1ub vis_seed_density;
		slimage::Image3ub vis_density_delta;
		if(plot_density_) {
			slimage::Image1f density = ComputeDepthDensity(clustering_.points, clustering_.opt);
			vis_density = slimage::Convert_f_2_ub(density, 20.0f);
			slimage::Image1f seed_density = ComputeDepthDensityFromSeeds(old_seeds, density, clustering_.opt);
			vis_seed_density = slimage::Convert_f_2_ub(seed_density, 20.0f);
			vis_density_delta.resize(density.width(), density.height());
			for(unsigned int i=0; i<density.getPixelCount(); i++) {
				vis_density_delta(i) = plots::PlusMinusColor(density(i) - seed_density(i), 0.025f);
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
		}

	}
	DANVIL_BENCHMARK_STOP(plotting)
}

void DaspTracker::trainInitialColorModel()
{
	const unsigned int cWidth = 640;
	const unsigned int cHeight = 480;

	if(boost::accumulators::count(coverage) > 0) {
		std::ofstream fs("coverage.txt", std::ios_base::app);
		double mean = boost::accumulators::mean(coverage);
		double var = boost::accumulators::variance(coverage);

		fs << mean << "\t" << std::sqrt(var) << std::endl;
		std::cout << "Mean=" << mean << ", sqrt(variance)=" << std::sqrt(var) << std::endl;
	}
	coverage = CoverageAccType();

	// capture all super pixels which are in the ROI and near to the camera
	LOG_NOTICE << "Initializing Appearance model --- STARTED";

	slimage::Image3ub initial_(cWidth, cHeight);
	initial_.fill(128);

	DANVIL_BENCHMARK_START(Hand_V3_Init)

	// find all clusters where the center is inside the roi
	const unsigned int R = 80;
	const unsigned int cRoiDepthBins = 100;
	const uint16_t cRoiDepthMin = 400;
	const uint16_t cRoiDepthMax = 2400;
	dasp::AutoFindDepthRange::DepthHistogram depth_hist(cRoiDepthMin, cRoiDepthMax, cRoiDepthBins);

	int roi_x1 = int(cWidth)/2 - R;
	int roi_x2 = int(cWidth)/2 + R;
	int roi_y1 = int(cHeight)/2 - R;
	int roi_y2 = int(cHeight)/2 + R;

	std::vector<dasp::Cluster> clusters_in_roi;
	std::vector<unsigned int> clusters_in_roi_ids;

	for(unsigned int k=0; k<clustering_.cluster.size(); k++) {
		const dasp::Cluster& cluster = clustering_.cluster[k];
		unsigned int in_roi_cnt = 0;
		for(unsigned int i : cluster.pixel_ids) {
			const dasp::Point& p = clustering_.points[i];
			if(roi_x1 <= p.spatial_x() && p.spatial_x() <= roi_x2
				&& roi_y1 <= p.spatial_y() && p.spatial_y() <= roi_y2) {
				in_roi_cnt ++;
			}
		}
		if(in_roi_cnt > cluster.pixel_ids.size()/2) {
			for(unsigned int i : cluster.pixel_ids) {
				const dasp::Point& p = clustering_.points[i];
				depth_hist.add(p.depth_i16);
			}
			clusters_in_roi.push_back(cluster);
			clusters_in_roi_ids.push_back(k);
		}
	}

	// find cut-off
	uint16_t depth_min;
	uint16_t depth_max;
	dasp::AutoFindDepthRange::Find(depth_hist, depth_min, depth_max);

	std::vector<dasp::Cluster> clusters_selected;
	std::vector<unsigned int> clusters_selected_id;
	for(unsigned int i=0; i<clusters_in_roi.size(); i++) {
		const dasp::Cluster& cluster = clusters_in_roi[i];
		uint16_t d = cluster.center.depth_i16;
		if(depth_min <= d && d <= depth_max) {
			clusters_selected.push_back(cluster);
			clusters_selected_id.push_back(clusters_in_roi_ids[i]);
		}
	}

	{
		// render all super pixel in range
		for(const dasp::Cluster& cluster : clusters_selected) {
			plots::PlotClusterPoints(initial_, cluster, clustering_.points, plots::RgbColor(cluster.center));
		}

//		// render roi
//		Danvil::Images::ImagePaint::PaintLine(initial_, roi_x1, roi_y1, roi_x2, roi_y1, Danvil::Images::Pixel3ub(255,0,0));
//		Danvil::Images::ImagePaint::PaintLine(initial_, roi_x1, roi_y2, roi_x2, roi_y2, Danvil::Images::Pixel3ub(255,0,0));
//		Danvil::Images::ImagePaint::PaintLine(initial_, roi_x1, roi_y1, roi_x1, roi_y2, Danvil::Images::Pixel3ub(255,0,0));
//		Danvil::Images::ImagePaint::PaintLine(initial_, roi_x2, roi_y1, roi_x2, roi_y2, Danvil::Images::Pixel3ub(255,0,0));
	}

	// learn gmm from contents of selected super pixels
	if(clusters_selected.size() >= 5) {
		std::vector<Danvil::ctLinAlg::Vec3f> gmm_training;
//		Superpixels6D::ForPixelClusters(clusters_selected, points,
//				[&gmm_training](unsigned int, const Superpixels6D::Cluster&, unsigned int, const Superpixels6D::Point& p) {
//			gmm_training.push_back(Danvil::ctLinAlg::Vec3f(p.color[0], p.color[1], p.color[2]));
//		});
		for(const dasp::Cluster& cluster : clusters_selected) {
//			gmm_training.push_back(Danvil::ctLinAlg::Vec3f(cluster.center.color[0], cluster.center.color[1], cluster.center.color[2]));
			for(unsigned int i : cluster.pixel_ids) {
				const dasp::Point& p = clustering_.points[i];
				gmm_training.push_back(Danvil::ctLinAlg::Vec3f(p.color[0], p.color[1], p.color[2]));
			}
		}
		hand_gmm_model_ = Danvil::GMM::GmmExpectationMaximization<2,3,float>(gmm_training);
		std::cout << hand_gmm_model_ << std::endl;
		std::cout << std::sqrt(std::abs(Danvil::ctLinAlg::Det(hand_gmm_model_.gaussians_[0].sigma()))) << std::endl;
		has_hand_gmm_model_ = true;

		// collect histograms for superpixels
		model_.reset(new SuperpixelHistogramModel());
		SuperpixelGraph G = clustering_.CreateNeighborhoodGraph();
		model_->train(G, clusters_selected_id);
	}

	DANVIL_BENCHMARK_STOP(Hand_V3_Init)

	{
		boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);

		images_["train"] = slimage::Ptr(initial_);
	}

	LOG_NOTICE << "Initializing Appearance model --- FINISHED";
}

void DaspTracker::performTrackingStep()
{

}

std::map<std::string, slimage::ImagePtr> DaspTracker::getImages() const
{
	boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);
	return images_;
}

slimage::Image1ub DaspTracker::getResultImage() const
{
	boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);
	return result_;
}

void DaspTracker::Render() const
{
	{	boost::interprocess::scoped_lock<boost::mutex> lock(render_mutex_);
		plots::RenderClusters(clustering_, cluster_color_mode_, selection_);
	}
}

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
