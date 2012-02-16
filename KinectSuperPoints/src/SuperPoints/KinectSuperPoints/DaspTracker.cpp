/*
 * DaspTracker.cpp
 *
 *  Created on: Feb 14, 2012
 *      Author: david
 */

#include "DaspTracker.h"
#include <SuperPoints/Mipmaps.hpp>
#include <SuperPoints/BlueNoise.hpp>
#include <SuperPoints/PointsAndNormals.hpp>
#include <SuperPoints/AutoDepth.hpp>
#define DANVIL_ENABLE_BENCHMARK
#include <Danvil/Tools/Benchmark.h>
#include <Danvil/Color.h>
#include <Danvil/Color/LAB.h>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <stdexcept>
using namespace std;
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

DaspTracker::DaspTracker()
{
	training_ = false;
	dasp_params.reset(new dasp::Parameters());
	dasp_params->focal = 580.0f;
	dasp_params->seed_mode = dasp::SeedModes::DepthMipmap;
	color_model_sigma_scale_ = 1.0f;
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

void DaspTracker::step(Danvil::Images::Image1ui16Ptr raw_kinect_depth, Danvil::Images::Image3ubPtr raw_kinect_color)
{
	DANVIL_BENCHMARK_START(step)

	// kinect 16-bit depth image
	kinect_depth.resize(raw_kinect_depth->width(), raw_kinect_depth->height());
	for(unsigned int i=0; i<kinect_depth.size(); i++) {
		uint16_t d = (*raw_kinect_depth)[i];
		if(d > 5000) {
			d = 0;
		}
		kinect_depth[i] = d;
	}

	// kinect RGB color image
	kinect_color_rgb.resize(raw_kinect_color->width(), raw_kinect_color->height());
	kinect_color_rgb.copyFrom(raw_kinect_color->begin());

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
		lab[0] = r;
		lab[1] = g;
		lab[2] = b;
	}, slimage::ThreadingOptions(1));

	DANVIL_BENCHMARK_STOP(step)

	performSegmentationStep();

	if(training_) {
		trainInitialColorModel();
		training_ = false;
	}

	performTrackingStep();

	DANVIL_BENCHMARK_PRINTALL_COUT
}

slimage::Image3ub ColorizeIntensity(const slimage::Image1f& I, float min, float max)
{
	slimage::Image3ub col(I.width(), I.height());
	if(I) {
		Danvil::ContinuousIntervalColorMapping<unsigned char, float> cm
				= Danvil::ContinuousIntervalColorMapping<unsigned char, float>::Factor_Blue_Red_Yellow_White();
		cm.useCustomBorderColors(Danvil::Colorub::Black, Danvil::Colorub::White);
		cm.setRange(min, max);
		slimage::ParallelProcess(I, col, [&cm](const float* src, unsigned char* dst) {
			cm(*src).writeRgb(dst);
		});
	}
	return col;
}

void DaspTracker::performSegmentationStep()
{
	// superpixel parameters
	dasp::ParametersExt super_params_ext = dasp::ComputeParameters(*dasp_params, kinect_color.width(), kinect_color.height());

	// compute normals only if necessary
	slimage::Image3f kinect_normals;
	if(super_params_ext.weight_normal > 0.0f) {
		slimage::Image3f kinect_points = dasp::PointsAndNormals::ComputePoints(kinect_depth);

		DANVIL_BENCHMARK_START(normals)
		slimage::Image3f kinect_normals = dasp::PointsAndNormals::ComputeNormals(kinect_depth, kinect_points);
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


		std::vector<float> cluster_probability_2 = dasp::ClassifyClusters(clusters,
				[&clusters, this](const dasp::Point& center) {
					dasp::SuperpixelHistogram hist(dasp::SuperpixelState({center.color, center.normal}));
					// find in range
					for(const dasp::Cluster& n : clusters) {
						float dmax = center.scala;
						if((center.pos - n.center.pos).norm() < dmax) {
							dasp::SuperpixelState n_state{n.center.color, Eigen::Vector3f::UnitZ()}; // TODO n.center.normal;
							hist.add(n_state);
						}
					}
					// compute minimum to learned clusters
					float c_min = 1e9;
					for(const dasp::SuperpixelHistogram& h : model_hist_) {
						float c = dasp::SuperpixelHistogram::Distance(hist, h);
						c_min = std::min(c, c_min);
					}
					std::cout << c_min << std::endl;
					return 1.0f - 0.1f *c_min;
				});

		probability_2.resize(kinect_depth.width(), kinect_depth.height());
		probability_2.fill(-1.0f);
		dasp::ForPixelClusters(clusters, points, [&probability_2,&cluster_probability_2](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
			float p_base = cluster_probability_2[cid];
			probability_2(p.spatial_x(), p.spatial_y()) = p_base;
		});

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
			});
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
			probability_color = ColorizeIntensity(probability, 0.0f, 1.0f);
			dasp::PlotEdges(superpixel_labels, probability_color, 2, 0, 0, 0);
		}

		if(probability_2) {
			probability_2_color = ColorizeIntensity(probability_2, 0.0f, 1.0f);
			dasp::PlotEdges(superpixel_labels, probability_2_color, 2, 0, 0, 0);
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

	for(const dasp::Cluster& cluster : clusters) {
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
		}
	}

	// find cut-off
	uint16_t depth_min;
	uint16_t depth_max;
	dasp::AutoFindDepthRange::Find(depth_hist, depth_min, depth_max);

	std::vector<dasp::Cluster> clusters_selected;
	for(const dasp::Cluster& cluster : clusters_in_roi) {
		uint16_t d = uint16_t(cluster.center.depth * 1000.0f); // TODO ........ pathetic
		if(depth_min <= d && d <= depth_max) {
			clusters_selected.push_back(cluster);
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
	if(clusters_selected.size() > 0) {
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
		model_hist_.clear();
		for(const dasp::Cluster& cluster : clusters_selected) {
			dasp::SuperpixelState state{cluster.center.color, Eigen::Vector3f::UnitZ()}; // TODO cluster.center.normal;
			dasp::SuperpixelHistogram hist(state);
			// iterate over neighbours
			for(const dasp::Cluster& n : clusters_selected) {
				//float dmax = std::max(cluster.center.scala, n.center.scala);
				float dmax = cluster.center.scala;
				if((cluster.center.pos - n.center.pos).norm() < dmax) {
					dasp::SuperpixelState n_state{n.center.color, Eigen::Vector3f::UnitZ()}; // TODO n.center.normal;
					hist.add(n_state);
				}
			}
			std::cout << "Chroma   : " << hist.hist_chroma_.transpose() << std::endl;
			std::cout << "Intensity: " << hist.hist_intensity_.transpose() << std::endl;
			std::cout << "Normal   : " << hist.hist_normal_.transpose() << std::endl;
			model_hist_.push_back(hist);
		}
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

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
