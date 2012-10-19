#include "WdgtKinectSuperPoints.h"
#include <dasp/IO.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Parallel.h>
#include <boost/bind.hpp>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <boost/format.hpp>

WdgtKinectSuperPoints::WdgtKinectSuperPoints(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);

	mode_ = IdleMode;
	capture_next_ = false;
	frame_counter_ = 0;
	has_new_frame_ = false;
	save_dasp_enabled_ = false;
	save_dasp_fn_ = "/tmp/dasp";

	QObject::connect(ui.pushButtonLoadOne, SIGNAL(clicked()), this, SLOT(OnLoadOne()));
#if defined DASP_HAS_OPENNI
	QObject::connect(ui.pushButtonSaveOne, SIGNAL(clicked()), this, SLOT(OnCaptureOne()));
	QObject::connect(ui.pushButtonLoadOni, SIGNAL(clicked()), this, SLOT(OnLoadOni()));
	QObject::connect(ui.pushButtonLive, SIGNAL(clicked()), this, SLOT(OnLive()));
#else
	ui.pushButtonSaveOne->setEnabled(false);
	ui.pushButtonLoadOni->setEnabled(false);
	ui.pushButtonLive->setEnabled(false);
#endif
	QObject::connect(ui.pushButtonSaveDebug, SIGNAL(clicked()), this, SLOT(OnSaveDebugImages()));
	QObject::connect(ui.pushButtonSaveDasp, SIGNAL(clicked()), this, SLOT(OnSaveDasp()));

	dasp_processing_.reset(new DaspProcessing());

	gui_params_.reset(new WdgtSuperpixelParameters(dasp_processing_));
	gui_params_->reload = &reload;
	gui_params_->show();

	gui_benchmark_.reset(new WdgtBenchmark());
	Danvil::Benchmark::Instance().setOnUpdate(boost::bind(&WdgtBenchmark::update, gui_benchmark_, _1, _2));
	ui.verticalLayoutBenchmark->addWidget(gui_benchmark_.get());

#if defined DASP_HAS_CANDY
	std::cout << "Creating OpenGL Widget ..." << std::endl;

	view_ = Candy::View::FactorDefaultPerspectiveView();
	scene_ = view_->getScene();

	boost::shared_ptr<Candy::DirectionalLight> light1(new Candy::DirectionalLight(Danvil::ctLinAlg::Vec3f(+1.0f, +1.0f, -1.0f)));
	light1->setDiffuse(Danvil::Colorf(1.0f, 1.0f, 1.0f));
	scene_->addLight(light1);
	boost::shared_ptr<Candy::DirectionalLight> light2(new Candy::DirectionalLight(Danvil::ctLinAlg::Vec3f(-1.0f, -1.0f, -1.0f)));
	light2->setDiffuse(Danvil::Colorf(1.0f, 1.0f, 1.0f));
	scene_->addLight(light2);

	engine_.reset(new Candy::Engine(view_));
	engine_->setClearColor(Danvil::Color::Grey);

	scene_->setShowCoordinateCross(false);

	gl_wdgt_ = new Candy::GLSystemQtWindow(0, engine_);
	ui.tabs->addTab(gl_wdgt_, "3D");

	boost::shared_ptr<Candy::IRenderable> dasp_renderling(
			new Candy::ObjectRenderling([this](){ dasp_processing_->Render(); }));
	scene_->addItem(dasp_renderling);
#endif

	QObject::connect(&timer_, SIGNAL(timeout()), this, SLOT(OnUpdateImages()));
	timer_.setInterval(1);
	timer_.start();

}

WdgtKinectSuperPoints::~WdgtKinectSuperPoints()
{
	StopProcessingThread();
}

void WdgtKinectSuperPoints::StopProcessingThread()
{
	// stops the processing thread and the kinect if it is running
#if defined DASP_HAS_OPENNI
	if(kinect_grabber_) {
		kinect_grabber_->Stop();
	}
#endif
	interrupt_loaded_thread_ = true;
	processing_thread_.join();
#if defined DASP_HAS_OPENNI
	kinect_grabber_.reset();
#endif
}

#if defined DASP_HAS_OPENNI

void WdgtKinectSuperPoints::OnCaptureOne()
{
	// this is only possible with a running kinect (ONI or live)
	if(!kinect_grabber_) {
		QMessageBox::information(this, "Not running!", "Open ONI file or run kinect in live mode!");
		return;
	}

	// user shall give rgbd image base name to save captured RGBD image
	QString fn = QFileDialog::getSaveFileName(this, "Open ONI file", "/home/david");
	if(fn.isEmpty()) {
		return;
	}
	capture_filename_ = fn.toStdString();

	// instruct to save the next image
	capture_next_ = true;
}

#endif

void WdgtKinectSuperPoints::OnLoadOne()
{
	// user shall select rgbd image base name
	// NOTE: Images are required to have the form BASE_color.png and BASE_depth.pgm,
	// NOTE: where BASE is the path and string selected in the dialog.
	QString fn = QFileDialog::getSaveFileName(this, "Open Color/Depth files", "/home/david");
	if(fn.isEmpty()) {
		return;
	}
	capture_filename_ = fn.toStdString();

	StopProcessingThread();

	// load selected images
	std::string fn_color = fn.toStdString() + "_color.png";
	std::string fn_depth = fn.toStdString() + "_depth.pgm";
	std::cout << "Reading '" << fn_color << "' and '" << fn_depth << "'..." << std::flush;
	slimage::Image3ub loaded_kinect_color = slimage::Load3ub(fn_color);
	slimage::Image1ui16 loaded_kinect_depth = slimage::Load1ui16(fn_depth);
	if(loaded_kinect_color.width() != loaded_kinect_depth.width() || loaded_kinect_color.height() != loaded_kinect_depth.height()) {
		std::cerr << "Size of color and depth image must match!" << std::endl;
		return;
	}

	mode_ = SingleFileMode;

	// processing thread repeats loaded image (only updated if requested)
	processing_thread_ = boost::thread([this,loaded_kinect_color,loaded_kinect_depth]() {
		this->interrupt_loaded_thread_ = false;
		while(!this->interrupt_loaded_thread_) {
			slimage::Image3ub cpy_color = loaded_kinect_color.clone();
			slimage::Image1ui16 cpy_depth = loaded_kinect_depth.clone();
			this->OnImages(cpy_depth, cpy_color);
			frame_counter_ = 1;
			while(!this->reload && !this->interrupt_loaded_thread_) {
				usleep(1000);
			}
			this->reload = false;
		}
	});

}

#if defined DASP_HAS_OPENNI

void WdgtKinectSuperPoints::OnLoadOni()
{
	// the user shall select an ONI file
	QString fn = QFileDialog::getOpenFileName(this, "Open ONI file", "/home/david", tr("ONI files (*.oni)"));
	if(fn.isEmpty()) {
		return;
	}
	LoadOni(fn.toStdString());
}

void WdgtKinectSuperPoints::LoadOni(const std::string& fn)
{
	StopProcessingThread();

	// open Kinect OpenNI ONI file
	std::cout << "Opening oni file: " << fn << std::endl;
	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->OpenFile(fn);
	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));

	mode_ = ReplayOniMode;
	remembered_capture_fn_ = fn;
	
	// processing thread polls kinect
	frame_counter_ = 0;
	processing_thread_ = boost::thread(&Romeo::Kinect::KinectGrabber::Run, kinect_grabber_);
}

void WdgtKinectSuperPoints::OnLive()
{
	StopProcessingThread();

	// open Kinect in live mode
	QString config = "/home/david/Programs/RGBD/OpenNI/Platform/Linux-x86/Redist/Samples/Config/SamplesConfig.xml";
	std::cout << "Opening kinect config file: " << config.toStdString() << std::endl;
	kinect_grabber_.reset(new Romeo::Kinect::KinectGrabber());
	kinect_grabber_->OpenConfig(config.toStdString());
	kinect_grabber_->on_depth_and_color_.connect(boost::bind(&WdgtKinectSuperPoints::OnImages, this, _1, _2));

	mode_ = LiveMode;

	// processing thread polls kinect
	frame_counter_ = 0;
	processing_thread_ = boost::thread(&Romeo::Kinect::KinectGrabber::Run, kinect_grabber_);
}

#endif

void WdgtKinectSuperPoints::OnSaveDebugImages()
{
	// save all generated dasp images
	for(auto p : dasp_processing_->getImages()) {
		slimage::Save(p.second, "/tmp/dasp_" + p.first + ".png");
	}
}

void WdgtKinectSuperPoints::OnSaveDasp()
{
	if(ui.pushButtonSaveDasp->isChecked()) {
		// get filename
		QString fn = QFileDialog::getSaveFileName(this,
				"Choose filename base for saving images",
				QString::fromStdString(save_dasp_fn_));
		save_dasp_fn_ = fn.toStdString();
		save_dasp_enabled_ = true;
		this->reload = true; // trigger OnImages if in single file mode
		if(mode_ == ReplayOniMode) {
			// reload oni file because we want to record from the beginning
			// TODO make this optional
			// TODO use seeking
			LoadOni(remembered_capture_fn_);
		}
	}
	else {
		save_dasp_enabled_ = false;
	}
}

void WdgtKinectSuperPoints::OnImages(const slimage::Image1ui16& kinect_depth, const slimage::Image3ub& kinect_color)
{
	// save RGBD image if requested
	if(capture_next_) {
		slimage::Save(kinect_color, capture_filename_ + "_color.png");
		slimage::Save(kinect_depth, capture_filename_ + "_depth.pgm");
		capture_next_ = false;
	}
	// process image
	dasp::SetRandomNumberSeed(0);
	dasp_processing_->step(kinect_depth, kinect_color);
	frame_counter_ ++;
	has_new_frame_ = true;
	// save
	if(save_dasp_enabled_) {
		std::string fn1 = (boost::format(save_dasp_fn_+"%1$05d.tsv") % frame_counter_).str();
		dasp_processing_->clustering_.SaveToFile(fn1, false);
		if(boost::num_vertices(dasp_processing_->Gnb) > 0) {
			std::string fn2 = (boost::format(save_dasp_fn_+"%1$05d_graph.txt") % frame_counter_).str();
			dasp::SaveGraph(dasp_processing_->Gnb, fn2);
		}
		else {
			if(frame_counter_ == 1) {
				std::cerr << "Graph not created -> not saving graph." << std::endl;
			}
		}
	}
}

void WdgtKinectSuperPoints::OnUpdateImages()
{
	if(!has_new_frame_) {
		return;
	}
	ui.labelFrameCounter->setText(QString("Frame: %1").arg(frame_counter_));
	std::map<std::string, bool> tabs_usage;
	std::map<std::string, QWidget*> tabs_labels;
	// prepare tabs
	for(int i=0; i<ui.tabs->count(); i++) {
		std::string text = ui.tabs->tabText(i).toStdString();
		tabs_usage[text] = (text == "3D");
		tabs_labels[text] = ui.tabs->widget(i);
	}
	// add tabs if necessary and display generated dasp images
	for(auto p : dasp_processing_->getImages()) {
		// image to display
		slimage::ImagePtr ref_img = p.second;
		if(!ref_img) {
			continue;
		}
		// convert to Qt image
		QImage* qimg = slimage::qt::ConvertToQt(ref_img);
		if(qimg == 0) {
			continue;
		}
		// find label and create new if none exists
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
		// display image in label
		qlabel->setPixmap(QPixmap::fromImage(*qimg));
		// cleanup
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
	has_new_frame_ = false;
}
