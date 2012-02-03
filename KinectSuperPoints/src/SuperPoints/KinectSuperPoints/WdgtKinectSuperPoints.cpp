#include "WdgtKinectSuperPoints.h"
#include <Slimage/Qt.hpp>

WdgtKinectSuperPoints::WdgtKinectSuperPoints(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);

	dasp_params.reset(new dasp::Parameters());

	gui_params_.reset(new WdgtSuperpixelParameters(dasp_params));
	gui_params_->show();

	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->options().EnableDepthRange(0.4, 2.4);

	//kinect_grabber_->OpenFile(settings.oni_file);
	kinect_grabber_->OpenConfig("/home/david/Programs/RGBD/OpenNI/Platform/Linux-x86/Redist/Samples/Config/SamplesConfig.xml");

	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));

	kinect_thread_ = boost::thread(&Romeo::Kinect::KinectGrabber::Run, kinect_grabber_);

	QObject::connect(&timer_, SIGNAL(timeout()), this, SLOT(OnUpdateImages()));
	timer_.setInterval(20);
	timer_.start();
}

WdgtKinectSuperPoints::~WdgtKinectSuperPoints()
{

}

void WdgtKinectSuperPoints::OnImages(Danvil::Images::Image1ui16Ptr raw_kinect_depth, Danvil::Images::Image3ubPtr raw_kinect_color)
{
	// copy images
	slimage::Image1ui16 kinect_depth(raw_kinect_depth->width(), raw_kinect_depth->height());
	kinect_depth.copyFrom(raw_kinect_depth->begin());

	slimage::Image3ub kinect_color(raw_kinect_color->width(), raw_kinect_color->height());
	kinect_color.copyFrom(raw_kinect_color->begin());

	images_mutex_.lock();
	images_["depth"] = slimage::Ptr(kinect_depth);
	images_["color"] = slimage::Ptr(kinect_color);
	images_mutex_.unlock();
}

void WdgtKinectSuperPoints::OnUpdateImages()
{
	std::map<std::string, bool> tabs_usage;
	std::map<std::string, QWidget*> tabs_labels;
	// prepare tabs
	for(int i=0; i<ui.tabs->count(); i++) {
		std::string text = ui.tabs->tabText(i).toStdString();
		tabs_usage[text] = false;
		tabs_labels[text] = ui.tabs->widget(i);
	}
	// add / renew tabs
	images_mutex_.lock();
	for(auto p : images_) {
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
		qimgscl = qimgscl.mirrored(false, true);
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
	images_mutex_.unlock();
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
