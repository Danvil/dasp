#include "WdgtKinectSuperPoints.h"
#include <Slimage/Qt.hpp>
#include <Slimage/impl/io.hpp>
#include <SuperPoints/Mipmaps.hpp>
#include <SuperPoints/BlueNoise.hpp>

WdgtKinectSuperPoints::WdgtKinectSuperPoints(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);

	dasp_params.reset(new dasp::Parameters());
	dasp_params->focal = 580.0f;
//	dasp_params->seed_mode = dasp::SeedModes::DepthDependentMipmap;
	dasp_params->seed_mode = dasp::SeedModes::BlueNoise;

	gui_params_.reset(new WdgtSuperpixelParameters(dasp_params));
	gui_params_->show();

	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->options().EnableDepthRange(0.4, 2.4);

	kinect_grabber_->OpenFile("/home/david/WualaDrive/Danvil/DataSets/2012-01-12 Kinect Hand Motions/01_UpDown_Move.oni");
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
	kinect_depth.copyFrom(raw_kinect_depth->begin());

	// kinect RGB color image
	slimage::Image3ub kinect_color(raw_kinect_color->width(), raw_kinect_color->height());
	kinect_color.copyFrom(raw_kinect_color->begin());

	// visualization of kinect depth image
	slimage::Image1ub kinect_depth_8(raw_kinect_depth->width(), raw_kinect_depth->height());
	for(size_t i=0; i<kinect_depth_8.size(); i++) {
		unsigned int d16 = kinect_depth[i];
		unsigned int d8 = std::min(d16 >> 4, 255u);
		kinect_depth_8[i] = (d16 == 0 ? 0 : 255 - d8);
//		std::cout << d16 << " " << d8 << " " << kinect_depth[i] << std::endl;
	}

	// superpixel parameters
	dasp::ParametersExt super_params_ext = dasp::ComputeParameters(*dasp_params, kinect_color.width(), kinect_color.height());

	slimage::Image3f kinect_normals;

	// prepare super pixel points
	dasp::ImagePoints points = dasp::CreatePoints(kinect_color, kinect_depth, kinect_normals, super_params_ext);

	// compute super pixel point edges
	slimage::Image1f edges;
	dasp::ComputeEdges(points, edges, super_params_ext);

	// compute super pixels
//	std::vector<dasp::Cluster> clusters = dasp::ComputeSuperpixels(points, edges, super_params_ext);
	std::vector<dasp::Seed> seeds = dasp::FindSeeds(points, super_params_ext);
	slimage::Image1ub seeds_img(points.width(), points.height());
	seeds_img.fill(255);
	for(dasp::Seed s : seeds) {
		seeds_img(s.x, s.y) = 0;
		if(1 <= s.x) {
			seeds_img(s.x-1, s.y) = 0;
		}
		if(s.x + 1 < seeds_img.width()) {
			seeds_img(s.x+1, s.y) = 0;
		}
		if(1 <= s.y) {
			seeds_img(s.x, s.y-1) = 0;
		}
		if(s.y + 1 < seeds_img.width()) {
			seeds_img(s.x, s.y+1) = 0;
		}
	}
	if(!edges.isNull()) {
		dasp::ImproveSeeds(seeds, points, edges, super_params_ext);
	}
	std::vector<dasp::Cluster> clusters = dasp::ComputeSuperpixels(points, seeds, super_params_ext);

	// super pixel visualization
	slimage::Image3ub super(points.width(), points.height());
	dasp::PlotCluster(clusters, points, super);
	std::vector<int> superpixel_labels = dasp::ComputePixelLabels(clusters, points);
	dasp::PlotEdges(superpixel_labels, super, 2, 0, 0, 0);

	slimage::Image1f num(points.width(), points.height());
	for(unsigned int i=0; i<points.size(); i++) {
		num[i] = points[i].estimatedCount();
	}
	std::vector<slimage::Image1f> mipmaps = dasp::Mipmaps::ComputeMipmaps(num, 16);

	// set images for gui
	images_mutex_.lock();
	images_["depth"] = slimage::Ptr(kinect_depth_8);
	images_["color"] = slimage::Ptr(kinect_color);
	images_["super"] = slimage::Ptr(super);
	images_["mm1"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[1], 256.0f));
	images_["mm2"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[2], 16.0f));
	images_["mm3"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[3], 4.0f));
	images_["mm4"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[4], 1.0f));
	images_["mm5"] = slimage::Ptr(slimage::Convert_f_2_ub(mipmaps[5], 0.25f));
	images_["seeds"] = slimage::Ptr(seeds_img);
	images_mutex_.unlock();
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
