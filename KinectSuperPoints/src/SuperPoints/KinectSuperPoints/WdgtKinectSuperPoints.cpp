#include "WdgtKinectSuperPoints.h"
#include <Slimage/Qt.hpp>
#include <Slimage/impl/io.hpp>
#include <Slimage/Parallel.h>

WdgtKinectSuperPoints::WdgtKinectSuperPoints(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);

	dasp_tracker_.reset(new dasp::DaspTracker());

	gui_params_.reset(new WdgtSuperpixelParameters(dasp_tracker_->dasp_params));
	gui_params_->on_train_ = [this]() { dasp_tracker_->training_ = true; };
	gui_params_->on_change_cm_sigma_scale_ = [this](float val) { dasp_tracker_->color_model_sigma_scale_ = val; };
	gui_params_->show();

	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->options().EnableDepthRange(0.4, 2.4);

//	kinect_grabber_->OpenFile("/home/david/WualaDrive/Danvil/DataSets/2012-01-12 Kinect Hand Motions/01_UpDown_Move.oni");
	kinect_grabber_->OpenConfig("/home/david/Programs/RGBD/OpenNI/Platform/Linux-x86/Redist/Samples/Config/SamplesConfig.xml");

	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));

	kinect_thread_ = boost::thread(&Romeo::Kinect::KinectGrabber::Run, kinect_grabber_);

	QObject::connect(&timer_, SIGNAL(timeout()), this, SLOT(OnUpdateImages()));
	timer_.setInterval(20);
	timer_.start();

}

WdgtKinectSuperPoints::~WdgtKinectSuperPoints()
{
	//kinect_grabber_->
	kinect_thread_.join();
}

void WdgtKinectSuperPoints::OnImages(Danvil::Images::Image1ui16Ptr raw_kinect_depth, Danvil::Images::Image3ubPtr raw_kinect_color)
{
	// kinect 16-bit depth image
	slimage::Image1ui16 kinect_depth;
	kinect_depth.resize(raw_kinect_depth->width(), raw_kinect_depth->height());
	for(unsigned int i=0; i<kinect_depth.size(); i++) {
		uint16_t d = (*raw_kinect_depth)[i];
		if(d > 5000) {
			d = 0;
		}
		kinect_depth[i] = d;
	}

	// kinect RGB color image
	slimage::Image3ub kinect_color_rgb;
	kinect_color_rgb.resize(raw_kinect_color->width(), raw_kinect_color->height());
	kinect_color_rgb.copyFrom(raw_kinect_color->begin());


	dasp_tracker_->step(kinect_depth, kinect_color_rgb);
}

void WdgtKinectSuperPoints::OnUpdateImages()
{
	images_mutex_.lock();
	std::map<std::string, slimage::ImagePtr> images_tmp = dasp_tracker_->getImages();
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
//		qlabel->setScaledContents(true);
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
