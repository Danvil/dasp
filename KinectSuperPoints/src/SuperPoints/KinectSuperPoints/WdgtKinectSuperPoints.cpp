#include "WdgtKinectSuperPoints.h"
#include <Slimage/Qt.hpp>
#include <Slimage/impl/io.hpp>
#include <Slimage/Parallel.h>
#include <SuperPoints/Mipmaps.hpp>
#include <SuperPoints/BlueNoise.hpp>
#include <SuperPoints/PointsAndNormals.hpp>
#define DANVIL_ENABLE_BENCHMARK
#include <Danvil/Tools/Benchmark.h>
#include <Danvil/Color.h>

WdgtKinectSuperPoints::WdgtKinectSuperPoints(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);

	dasp_params.reset(new dasp::Parameters());
	dasp_params->focal = 580.0f;
	dasp_params->seed_mode = dasp::SeedModes::DepthMipmap;

	gui_params_.reset(new WdgtSuperpixelParameters(dasp_params));
	gui_params_->show();

	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->options().EnableDepthRange(0.4, 2.4);

//	kinect_grabber_->OpenFile("/home/david/WualaDrive/Danvil/DataSets/2012-01-12 Kinect Hand Motions/01_UpDown_Move.oni");
	kinect_grabber_->OpenConfig("/home/david/Programs/RGBD/OpenNI/Platform/Linux-x86/Redist/Samples/Config/SamplesConfig.xml");

	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));

	kinect_thread_ = boost::thread(&Romeo::Kinect::KinectGrabber::Run, kinect_grabber_);

	running_ = true;
//	kinect_thread_ = boost::thread(&WdgtKinectSuperPoints::ComputeBlueNoiseImpl, this);

	QObject::connect(&timer_, SIGNAL(timeout()), this, SLOT(OnUpdateImages()));
	timer_.setInterval(20);
	timer_.start();
}

WdgtKinectSuperPoints::~WdgtKinectSuperPoints()
{
	running_ = false;
	kinect_thread_.join();
}

void WdgtKinectSuperPoints::ComputeBlueNoiseImpl()
{
	slimage::Image1ub img_raw = slimage::Load1ub("/home/david/Documents/DataSets/2012-02-06 Blue Noise/flower_256.png");
	slimage::Image1f density(img_raw.width(), img_raw.height());
	for(unsigned int i=0; i<density.size(); i++) {
		density[i] = 0.25f*(1.0f - float(img_raw[i]) / 255.0f);
	}

	while(running_) {
		std::vector<dasp::BlueNoise::Point> points = dasp::BlueNoise::Compute(density);
		slimage::Image1ub img_pnts(density.width(), density.height());
		img_pnts.fill(255);
		dasp::BlueNoise::PlotPoints(points, img_pnts);

		// set images for gui
		images_mutex_.lock();
		images_["raw"] = slimage::Ptr(img_raw);
		images_["blue"] = slimage::Ptr(img_pnts);
		images_mutex_.unlock();
	}
}

void WdgtKinectSuperPoints::OnImages(Danvil::Images::Image1ui16Ptr raw_kinect_depth, Danvil::Images::Image3ubPtr raw_kinect_color)
{
	// kinect 16-bit depth image
	slimage::Image1ui16 kinect_depth(raw_kinect_depth->width(), raw_kinect_depth->height());
	for(unsigned int i=0; i<kinect_depth.size(); i++) {
		uint16_t d = (*raw_kinect_depth)[i];
//		if(d > 2000) {
//			d = 0;
//		}
		kinect_depth[i] = d;
	}

	// kinect RGB color image
	slimage::Image3ub kinect_color(raw_kinect_color->width(), raw_kinect_color->height());
	kinect_color.copyFrom(raw_kinect_color->begin());

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
	dasp::ImagePoints points = dasp::CreatePoints(kinect_color, kinect_depth, kinect_normals, super_params_ext);
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
	std::vector<dasp::Cluster> clusters = dasp::ComputeSuperpixels(points, seeds, super_params_ext);
	DANVIL_BENCHMARK_STOP(clusters)

	{
		DANVIL_BENCHMARK_START(plotting)

		// visualization of kinect depth image
		slimage::Image3ub kinect_depth_color(raw_kinect_depth->width(), raw_kinect_depth->height());
		for(size_t i=0; i<kinect_depth_color.size()/3; i++) {
			unsigned int d16 = kinect_depth[i];
			slimage::Pixel3ub color;
			if(d16 == 0) {
				kinect_depth_color(i,0) = {{0,0,0}};
			}
			else {
				// blue -> red -> yellow
				auto cm = Danvil::ContinuousIntervalColorMapping<unsigned char, uint16_t>::Factor_Blue_Red_Yellow();
				cm.setRange(400,2000);
				auto color = cm(d16);
				unsigned int q = d16 % 25;
				color.r = std::max(0, int(color.r) - int(q));
				color.g = std::max(0, int(color.g) - int(q));
				color.b = std::max(0, int(color.b) - int(q));
				color.writeRgb(kinect_depth_color.pointer(i,0));
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
//		std::vector<int> superpixel_labels = dasp::ComputePixelLabels(clusters, points);
//		dasp::PlotEdges(superpixel_labels, super, 2, 0, 0, 0);

	//	// plot mipmaps
	//	slimage::Image1f num(points.width(), points.height());
	//	for(unsigned int i=0; i<points.size(); i++) {
	//		num[i] = points[i].estimatedCount();
	//	}
	//	std::vector<slimage::Image1f> mipmaps = dasp::Mipmaps::ComputeMipmaps(num, 16);

		// set images for gui
		images_mutex_.lock();
		images_["color"] = slimage::Ptr(kinect_color);
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
		images_mutex_.unlock();

		DANVIL_BENCHMARK_STOP(clusters)
	}

	DANVIL_BENCHMARK_PRINTALL_COUT
}

void WdgtKinectSuperPoints::OnUpdateImages()
{
	images_mutex_.lock();
	std::map<std::string, slimage::ImagePtr> images_tmp = images_;
	images_mutex_.unlock();

	std::map<std::string, bool> tabs_usage;
	std::map<std::string, QWidget*> tabs_labels;
	// prepare tabs
	for(int i=0; i<ui.tabs->count(); i++) {
		std::string text = ui.tabs->tabText(i).toStdString();
		tabs_usage[text] = false;
		tabs_labels[text] = ui.tabs->widget(i);
	}
	// add / renew tabs
	for(auto p : images_tmp) {
		slimage::ImagePtr ref_img = p.second;
		if(!ref_img) {
			continue;
		}

		QImage* qimg = slimage::ConvertToQt(ref_img);
		if(qimg == 0) {
			continue;
		}

		QImage qimgscl;
//		if(qimg->width() > cFeedbackVisualSize || qimg->height() > cFeedbackVisualSize) {
//			qimgscl = qimg->scaled(QSize(cFeedbackVisualSize,cFeedbackVisualSize), Qt::KeepAspectRatio);
//		}
//		else {
			qimgscl = *qimg;
//		}
//		qimgscl = qimgscl.mirrored(false, true);
		QLabel* qlabel;
		std::string text = p.first;
		if(tabs_labels.find(text) != tabs_labels.end()) {
			qlabel = (QLabel*)tabs_labels[text];
		} else {
			qlabel = new QLabel();
			tabs_labels[text] = qlabel;
			ui.tabs->addTab(qlabel, QString::fromStdString(text));
		}
		tabs_usage[text] = true;
		qlabel->setPixmap(QPixmap::fromImage(qimgscl));
		delete qimg;
	}
	// delete unused tabs
	for(auto p : tabs_usage) {
		if(!p.second) {
			for(int i=0; i<ui.tabs->count(); i++) {
				if(ui.tabs->tabText(i).toStdString() == p.first) {
					ui.tabs->removeTab(i);
				}
			}
		}
	}
}
