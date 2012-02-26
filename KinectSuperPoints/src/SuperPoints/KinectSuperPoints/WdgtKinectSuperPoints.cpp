#include "WdgtKinectSuperPoints.h"
#include <Slimage/Qt.hpp>
#include <Slimage/impl/io.hpp>
#include <Slimage/Parallel.h>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>

WdgtKinectSuperPoints::WdgtKinectSuperPoints(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);

	QObject::connect(ui.pushButtonSaveOne, SIGNAL(clicked()), this, SLOT(OnCaptureOne()));
	QObject::connect(ui.pushButtonLoadOne, SIGNAL(clicked()), this, SLOT(OnLoadOne()));
	QObject::connect(ui.pushButtonLoadOni, SIGNAL(clicked()), this, SLOT(OnLoadOni()));
	QObject::connect(ui.pushButtonLive, SIGNAL(clicked()), this, SLOT(OnLive()));

	dasp_tracker_.reset(new dasp::DaspTracker());

	gui_params_.reset(new WdgtSuperpixelParameters(dasp_tracker_));
	gui_params_->show();

	LOG_NOTICE << "Creating OpenGL Widget ...";

	view_ = Danvil::SimpleEngine::View::FactorDefaultPerspectiveView();
	scene_ = view_->getScene();

	boost::shared_ptr<Danvil::SimpleEngine::DirectionalLight> light1(new Danvil::SimpleEngine::DirectionalLight(Danvil::ctLinAlg::Vec3f(+1.0f, +1.0f, -1.0f)));
	light1->setDiffuse(Danvil::Colorf(1.0f, 1.0f, 1.0f));
	view_->getScene()->addLight(light1);
	boost::shared_ptr<Danvil::SimpleEngine::DirectionalLight> light2(new Danvil::SimpleEngine::DirectionalLight(Danvil::ctLinAlg::Vec3f(-1.0f, -1.0f, -1.0f)));
	light2->setDiffuse(Danvil::Colorf(1.0f, 1.0f, 1.0f));
	view_->getScene()->addLight(light2);

	engine_.reset(new Danvil::SimpleEngine::Engine(view_));
	engine_->setClearColor(Danvil::Color::Grey);

	gl_wdgt_ = new Danvil::SimpleEngine::GLSystemQtWindow(0, engine_);
	ui.tabs->addTab(gl_wdgt_, "3D");

	boost::shared_ptr<Danvil::SimpleEngine::IRenderable> dasp_renderling(
			new Danvil::SimpleEngine::ObjectRenderling([this](){ dasp_tracker_->Render(); }));
	scene_->addItem(dasp_renderling);

	LOG_NOTICE << "Starting Kinect ...";

//	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
//	kinect_grabber_->options().EnableDepthRange(0.4, 2.4);
//	kinect_grabber_->OpenFile("/home/david/WualaDrive/Danvil/DataSets/2012-01-12 Kinect Hand Motions/01_UpDown_Move.oni");
//	kinect_grabber_->OpenFile("/home/david/Desktop/Kinect/2012-02-25 Stuff On Table.oni");
//	kinect_grabber_->OpenConfig("/home/david/Programs/RGBD/OpenNI/Platform/Linux-x86/Redist/Samples/Config/SamplesConfig.xml");
//	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));

	QObject::connect(&timer_, SIGNAL(timeout()), this, SLOT(OnUpdateImages()));
	timer_.setInterval(20);
	timer_.start();

}

WdgtKinectSuperPoints::~WdgtKinectSuperPoints()
{
}

void WdgtKinectSuperPoints::OnCaptureOne()
{
	if(!kinect_grabber_) {
		QMessageBox::information(this, "Not running!", "Open ONI file or run kinect in live mode!");
		return;
	}
	// FIXME
}

void WdgtKinectSuperPoints::OnLoadOne()
{
	if(!kinect_grabber_) {
		QMessageBox::information(this, "Not running!", "Open ONI file or run kinect in live mode!");
		return;
	}
	// FIXME
}

void WdgtKinectSuperPoints::OnLoadOni()
{
	QString fn = QFileDialog::getOpenFileName(this, "Open ONI file", "/home/david", tr("ONI files (*.oni)"));
	if(fn.isEmpty()) {
		return;
	}

	if(kinect_grabber_) {
		kinect_grabber_->stop_requested_ = true;
		kinect_thread_.join();
		kinect_grabber_.reset();
	}

	LOG_NOTICE << "Opening oni file: " << fn.toStdString();
	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->options().EnableDepthRange(0.4, 2.4);
	kinect_grabber_->OpenFile(fn.toStdString());
	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));
	kinect_thread_ = boost::thread(&Romeo::Kinect::KinectGrabber::Run, kinect_grabber_);
}

void WdgtKinectSuperPoints::OnLive()
{
	if(kinect_grabber_) {
		kinect_grabber_->stop_requested_ = true;
		kinect_thread_.join();
		kinect_grabber_.reset();
	}

	QString config = "/home/david/Programs/RGBD/OpenNI/Platform/Linux-x86/Redist/Samples/Config/SamplesConfig.xml";
	LOG_NOTICE << "Opening kinect config file: " << config.toStdString();
	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->options().EnableDepthRange(0.4, 2.4);
	kinect_grabber_->OpenConfig(config.toStdString());
	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));
	kinect_thread_ = boost::thread(&Romeo::Kinect::KinectGrabber::Run, kinect_grabber_);
}

void WdgtKinectSuperPoints::OnImages(Danvil::Images::Image1ui16Ptr raw_kinect_depth, Danvil::Images::Image3ubPtr raw_kinect_color)
{
	// kinect 16-bit depth image
	slimage::Image1ui16 kinect_depth;
	kinect_depth.resize(raw_kinect_depth->width(), raw_kinect_depth->height());
	for(unsigned int i=0; i<kinect_depth.size(); i++) {
		uint16_t d = (*raw_kinect_depth)[i];
		if(d > 3000) {
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
		tabs_usage[text] = (text == "3D");
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
