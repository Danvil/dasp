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
#define DANVIL_ENABLE_BENCHMARK
#include <Slimage/Paint.hpp>
#include <Danvil/Tools/Benchmark.h>
#include <Danvil/Color.h>
#include <Danvil/Color/LAB.h>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <stdexcept>
using namespace std;
//----------------------------------------------------------------------------//
namespace dasp {
//----------------------------------------------------------------------------//

DaspTracker::DaspTracker()
{
	training_ = false;
	dasp_params.reset(new dasp::Parameters());
	dasp_params->focal = 580.0f;
	dasp_params->seed_mode = dasp::SeedModes::DepthMipmap;
	color_model_sigma_scale_ = 1.0f;
	thread_pool_index_ = 100;
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

	kinect_depth = raw_kinect_depth.clone();

	kinect_color_rgb = raw_kinect_color.clone();

	// convert rgb to lab
	kinect_color.resize(kinect_color_rgb.width(), kinect_color_rgb.height());
	slimage::ParallelProcess(kinect_color_rgb, kinect_color, [](const unsigned char* rgb, float* lab) {
		float r = float(rgb[0]) / 255.0f;
		float g = float(rgb[1]) / 255.0f;
		float b = float(rgb[2]) / 255.0f;
//		float lab1, lab2, lab3;
//		Danvil::color_rgb_to_lab(r, g, b, lab1, lab2, lab3);
//		lab1 /= 100.0f;
//		lab2 /= 100.0f;
//		lab3 /= 100.0f;
////		std::cout << r << " " << g << " " << b << " " << lab1 << " " << lab2 << " " << lab3 << std::endl;
//		lab[0] = lab1;
//		lab[1] = lab2;
//		lab[2] = lab3;
		float a = r + g + b;
		if(a > 0.05f) {
			r /= a;
			g /= a;
			b /= a;
		}
		else {
			r = 0;
			g = 0;
			b = 0;
		}
		lab[0] = r;
		lab[1] = g;
		lab[2] = b;
	}, slimage::ThreadingOptions::UsePool(thread_pool_index_));

	DANVIL_BENCHMARK_STOP(step)

	performSegmentationStep();

	if(training_) {
		trainInitialColorModel();
		training_ = false;
	}

	performTrackingStep();

	DANVIL_BENCHMARK_PRINTALL_COUT
}

slimage::Image3ub ColorizeIntensity(const slimage::Image1f& I, float min, float max, unsigned int pool_id)
{
	slimage::Image3ub col(I.width(), I.height());
	if(I) {
		Danvil::ContinuousIntervalColorMapping<unsigned char, float> cm
				= Danvil::ContinuousIntervalColorMapping<unsigned char, float>::Factor_Blue_Red_Yellow_White();
		cm.useCustomBorderColors(Danvil::Colorub::Black, Danvil::Colorub::White);
		cm.setRange(min, max);
		slimage::ParallelProcess(I, col, [&cm](const float* src, unsigned char* dst) {
			cm(*src).writeRgb(dst);
		}, slimage::ThreadingOptions::UsePool(pool_id));
	}
	return col;
}

SuperpixelGraph CreateGraphFromClusters(const std::vector<dasp::Cluster>& clusters)
{
	SuperpixelGraph G;
	for(const dasp::Cluster& c : clusters) {
		SuperpixelState s;
		s.x = c.center.spatial_x();
		s.y = c.center.spatial_y();
		s.color = c.center.color;
		s.normal = c.center.normal;
		s.position = PointsAndNormals::PointFromZ(
				c.center.spatial_x(), c.center.spatial_y(),
				c.center.depth,
				640, 480);
		s.scala = c.center.scala;
		G.nodes_.push_back(s);
	}
	G.createConnections(0.10f); // FIXME constant !!!
	return G;
}

void DaspTracker::performSegmentationStep()
{
	// superpixel parameters
	dasp::ParametersExt super_params_ext = dasp::ComputeParameters(*dasp_params, kinect_color.width(), kinect_color.height());

	// compute normals only if necessary
	slimage::Image3f kinect_normals;
	if(super_params_ext.weight_normal > 0.0f) {
		slimage::Image3f kinect_points = dasp::PointsAndNormals::ComputePoints(kinect_depth, slimage::ThreadingOptions::UsePool(thread_pool_index_));

		DANVIL_BENCHMARK_START(normals)
		slimage::Image3f kinect_normals = dasp::PointsAndNormals::ComputeNormals(kinect_depth, kinect_points, slimage::ThreadingOptions::UsePool(thread_pool_index_));
	//	dasp::PointsAndNormals::ComputeNormalsFast(kinect_depth, kinect_points, kinect_normals);
		DANVIL_BENCHMARK_STOP(normals)
	}

	// prepare super pixel points
	DANVIL_BENCHMARK_START(points)
	points = dasp::CreatePoints(kinect_color, kinect_depth, kinect_normals, super_params_ext);
	DANVIL_BENCHMARK_STOP(points)

	// compute super pixel seeds
	DANVIL_BENCHMARK_START(seeds)
	std::vector<dasp::Seed> seeds = dasp::FindSeeds(points, super_params_ext);
	std::cout << "Seeds: " << seeds.size() << std::endl;
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
	clusters = dasp::ComputeSuperpixels(points, seeds, super_params_ext);
	DANVIL_BENCHMARK_STOP(clusters)

	slimage::Image1f probability;
	slimage::Image1f probability_2;
	slimage::Image3ub plot_graph;
	slimage::Image3ub plot_labels;

	if(has_hand_gmm_model_) {

		auto scaled_model = hand_gmm_model_;
		scaled_model.ScaleDeviation(color_model_sigma_scale_); // FIXME get parameter from Gui

		// classify color with gmm
		std::vector<float> cluster_probability = dasp::ClassifyClusters(clusters,
				[&scaled_model](const dasp::Point& center) {
					Danvil::ctLinAlg::Vec3f v(center.color[0], center.color[1], center.color[2]);
					return scaled_model(v);
				});

		// paint cluster probabilities
		probability.resize(kinect_depth.width(), kinect_depth.height());
		probability.fill(-1.0f);
		dasp::ForPixelClusters(clusters, points, [&probability,&cluster_probability,&scaled_model](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
//			Danvil::ctLinAlg::Vec3f v(p.color[0], p.color[1], p.color[2]);
//			float p_pixel = scaled_model(v);
			float p_base = cluster_probability[cid];
//			float p_final = 0.5f * (p_pixel + p_base);
			probability(p.spatial_x(), p.spatial_y()) = p_base;
		});

		// histogram color model
		SuperpixelGraph G = CreateGraphFromClusters(clusters);
		std::vector<float> cluster_probability_2 = model_->evaluate(G);

		probability_2.resize(kinect_depth.width(), kinect_depth.height());
		probability_2.fill(-1.0f);
		dasp::ForPixelClusters(clusters, points, [&probability_2,&cluster_probability_2](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
			float p_base = cluster_probability_2[cid];
			probability_2(p.spatial_x(), p.spatial_y()) = p_base;
		});

		plot_graph.resize(kinect_color.width(), kinect_color.height());
		plot_graph.fill(0);
		for(unsigned int i=0; i<G.nodes_.size(); i++) {
			for(unsigned int k : G.node_connections_[i]) {
				slimage::PaintLine(plot_graph, G.nodes_[i].x, G.nodes_[i].y, G.nodes_[k].x, G.nodes_[k].y, slimage::Pixel3ub{255,255,255});
			}
		}

		// plot model gaussian cluster id for each superpixel
		std::vector<unsigned int> labels = model_->label(G);
		plot_labels.resize(kinect_color.width(), kinect_color.height());
		plot_labels.fill(0);
		std::vector<Eigen::Vector3f> model_colors = model_->getClusterColors();
		dasp::ForPixelClusters(clusters, points, [&plot_labels,&labels,&model_colors](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
			unsigned int l = labels[cid];
			Eigen::Vector3f color = model_colors[l];
			plot_labels(p.spatial_x(), p.spatial_y()) = slimage::Pixel3ub{
				{(unsigned char)(255.0f*color[0]), (unsigned char)(255.0f*color[1]), (unsigned char)(255.0f*color[2])}
			};
		});


	}

	result_.resize(kinect_color.width(), kinect_color.height());
	if(has_hand_gmm_model_) {
		result_ = slimage::Convert_f_2_ub(probability);
	}
	else {
		result_.fill(0);
	}

	{
		DANVIL_BENCHMARK_START(plotting)

		// visualization of kinect depth image
		slimage::Image3ub kinect_depth_color(kinect_depth.width(), kinect_depth.height());
		for(size_t i=0; i<kinect_depth_color.size()/3; i++) {
			unsigned int d16 = kinect_depth[i];
			slimage::Pixel3ub color;
			if(d16 == 0) {
				kinect_depth_color(i) = {{0,0,0}};
			}
			else {
				// blue -> red -> yellow
				auto cm = Danvil::ContinuousIntervalColorMapping<unsigned char, uint16_t>::Factor_Blue_Red_Yellow();
				cm.setRange(400,2000);
				auto color = cm(d16);
				unsigned int q = d16 % 25;
				unsigned char r = std::max(0, int(color.r) - int(q));
				unsigned char g = std::max(0, int(color.g) - int(q));
				unsigned char b = std::max(0, int(color.b) - int(q));
				kinect_depth_color(i) = {{r,g,b}};
			}
		}

		// plot normals
		slimage::Image3ub kinect_normals_vis;
		if(kinect_normals) {
			kinect_normals_vis.resize(kinect_normals.width(), kinect_normals.height());
			slimage::ParallelProcess(kinect_normals, kinect_normals_vis, [](const float* src, unsigned char* dst) {
				dst[0] = int(128.0f + 128.0f * src[0]);
				dst[1] = int(128.0f + 128.0f * src[1]);
				dst[2] = int(128.0f + 128.0f * src[2]);
			}, slimage::ThreadingOptions::UsePool(thread_pool_index_));
		}

		// plot seeds
		slimage::Image1ub seeds_img(points.width(), points.height());
		seeds_img.fill(255);
		dasp::PlotSeeds(seeds, seeds_img, 0, false);

		// plot super pixel with edges
		slimage::Image3ub super(points.width(), points.height());
		dasp::PlotCluster(clusters, points, super);
		std::vector<int> superpixel_labels = dasp::ComputePixelLabels(clusters, points);
		dasp::PlotEdges(superpixel_labels, super, 2, 0, 0, 0);

	//	// plot mipmaps
	//	slimage::Image1f num(points.width(), points.height());
	//	for(unsigned int i=0; i<points.size(); i++) {
	//		num[i] = points[i].estimatedCount();
	//	}
	//	std::vector<slimage::Image1f> mipmaps = dasp::Mipmaps::ComputeMipmaps(num, 16);

		slimage::Image3ub probability_color;
		slimage::Image3ub probability_2_color;

		if(probability) {
			probability_color = ColorizeIntensity(probability, 0.0f, 1.0f, thread_pool_index_);
			//dasp::PlotEdges(superpixel_labels, probability_color, 2, 0, 0, 0);

//			for(unsigned int i=0; i<probability_color.size(); i++) {
//				if(probability(i) > 0.90f) {
//					probability_color(i)[0] = kinect_color_rgb(i)[0];
//					probability_color(i)[1] = kinect_color_rgb(i)[1];
//					probability_color(i)[2] = kinect_color_rgb(i)[2];
//				}
//			}

		}

		if(probability_2) {
			probability_2_color = ColorizeIntensity(probability_2, 0.0f, 1.0f, thread_pool_index_);
			//dasp::PlotEdges(superpixel_labels, probability_2_color, 2, 0, 0, 0);
		}

		// set images for gui
		{
			boost::interprocess::scoped_lock<boost::mutex> lock(images_mutex_);

			images_["rgb"] = slimage::Ptr(kinect_color_rgb);
//			images_["lab"] = slimage::Ptr(slimage::Convert_f_2_ub(kinect_color));
			images_["depth"] = slimage::Ptr(kinect_depth_color);
			images_["super"] = slimage::Ptr(super);
		//	images_["mm1"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[1], 256.0f));
		//	images_["mm2"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[2], 16.0f));
		//	images_["mm3"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[3], 4.0f));
		//	images_["mm4"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[4], 1.0f));
		//	images_["mm5"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[5], 0.25f));
			images_["seeds"] = slimage::Ptr(seeds_img);
			if(kinect_normals_vis) {
				images_["normals"] = slimage::Ptr(kinect_normals_vis);
			}
			images_["prob"] = slimage::Ptr(probability_color);
			images_["prob2"] = slimage::Ptr(probability_2_color);
			images_["graph"] = slimage::Ptr(plot_graph);
			images_["labels"] = slimage::Ptr(plot_labels);
		}

		DANVIL_BENCHMARK_STOP(clusters)
	}
}

void DaspTracker::trainInitialColorModel()
{
	const unsigned int cWidth = 640;
	const unsigned int cHeight = 480;

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

	for(unsigned int k=0; k<clusters.size(); k++) {
		const dasp::Cluster& cluster = clusters[k];
		unsigned int in_roi_cnt = 0;
		for(unsigned int i : cluster.pixel_ids) {
			const dasp::Point& p = points[i];
			if(roi_x1 <= p.spatial_x() && p.spatial_x() <= roi_x2
				&& roi_y1 <= p.spatial_y() && p.spatial_y() <= roi_y2) {
				in_roi_cnt ++;
			}
		}
		if(in_roi_cnt > cluster.pixel_ids.size()/2) {
			for(unsigned int i : cluster.pixel_ids) {
				const dasp::Point& p = points[i];
				uint16_t d = uint16_t(p.depth * 1000.0f); // TODO ........ pathetic
				depth_hist.add(d);
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
		uint16_t d = uint16_t(cluster.center.depth * 1000.0f); // TODO ........ pathetic
		if(depth_min <= d && d <= depth_max) {
			clusters_selected.push_back(cluster);
			clusters_selected_id.push_back(clusters_in_roi_ids[i]);
		}
	}

	{
		// render all super pixel in range
		for(const dasp::Cluster& cluster : clusters_selected) {
			PlotCluster(cluster, points, initial_);
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
				const dasp::Point& p = points[i];
				gmm_training.push_back(Danvil::ctLinAlg::Vec3f(p.color[0], p.color[1], p.color[2]));
			}
		}
		hand_gmm_model_ = Danvil::GMM::GmmExpectationMaximization<5,3,float>(gmm_training);
		std::cout << hand_gmm_model_ << std::endl;
		std::cout << std::sqrt(std::abs(Danvil::ctLinAlg::Det(hand_gmm_model_.gaussians_[0].sigma()))) << std::endl;
		has_hand_gmm_model_ = true;

		// collect histograms for superpixels
		model_.reset(new SuperpixelHistogramModel());
		SuperpixelGraph G = CreateGraphFromClusters(clusters);
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

//----------------------------------------------------------------------------//
}
//----------------------------------------------------------------------------//
