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

	QObject::connect(ui.horizontalSliderFrame, SIGNAL(valueChanged(int)), this, SLOT(OnChangeFrame(int)));
	ui.horizontalSliderFrame->setEnabled(false);
	ui.horizontalSliderFrame->setMinimum(0);
	ui.horizontalSliderFrame->setMaximum(0);
	QObject::connect(ui.actionLoad_RGB_D_Image, SIGNAL(triggered()), this, SLOT(OnLoadOne()));
#if defined DASP_HAS_OPENNI
	QObject::connect(ui.actionKinect_Live_Mode, SIGNAL(triggered()), this, SLOT(OnLive()));
	QObject::connect(ui.actionLoad_ONI, SIGNAL(triggered()), this, SLOT(OnLoadOni()));
	QObject::connect(ui.actionSave_RBG_D_Image, SIGNAL(triggered()), this, SLOT(OnCaptureOne()));
#else
	ui.actionKinect_Live_Mode->setEnabled(false);
	ui.actionLoad_ONI->setEnabled(false);
	ui.actionSave_RBG_D_Image->setEnabled(false);
#endif
	QObject::connect(ui.actionSave_Debug_Images, SIGNAL(triggered()), this, SLOT(OnSaveDebugImages()));
	QObject::connect(ui.actionBatch_Save, SIGNAL(triggered()), this, SLOT(OnSaveDasp()));

	QObject::connect(ui.action_Parameters, SIGNAL(triggered()), this, SLOT(onViewParameters()));
	QObject::connect(ui.action_Benchmarks, SIGNAL(triggered()), this, SLOT(onViewBenchmark()));
	QObject::connect(ui.actionAbout, SIGNAL(triggered()), this, SLOT(onViewAbout()));

	dasp_processing_.reset(new DaspProcessing());

	gui_params_.reset(new WdgtSuperpixelParameters(dasp_processing_));
	gui_params_->setAttribute(Qt::WA_DeleteOnClose, false);
	gui_params_->reload = &reload;
	gui_params_->show();

	gui_benchmark_.reset(new WdgtBenchmark());
	gui_benchmark_->setAttribute(Qt::WA_DeleteOnClose, false);
	Danvil::Benchmark::Instance().setOnUpdate(boost::bind(&WdgtBenchmark::update, gui_benchmark_, _1, _2));

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

void WdgtKinectSuperPoints::ShowLive()
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
	
	ui.horizontalSliderFrame->setEnabled(true);
	ui.horizontalSliderFrame->setMaximum(kinect_grabber_->NumFrames());

	// processing thread polls kinect
	frame_counter_ = 0;
	processing_thread_ = boost::thread(
		[this]() {
			this->interrupt_loaded_thread_ = false;
			while(!this->interrupt_loaded_thread_) {
				unsigned int requested_frame = frame_counter_ + (this->ui.checkBoxPlay->isChecked() ? 1 : 0);
				while(!this->reload && requested_frame == frame_counter_ && !this->interrupt_loaded_thread_) {
					usleep(1000);
					requested_frame = ui.horizontalSliderFrame->value() + (this->ui.checkBoxPlay->isChecked() ? 1 : 0);
				}
				this->reload = false;
				// seek to requested frame (playing will start from there)
				kinect_grabber_->SeekToFrame(requested_frame);
				frame_counter_ = kinect_grabber_->TellFrame();
				// grab frame from oni
				bool ok = this->kinect_grabber_->Grab();
				if(!ok) {
					break;
				}
				// update images with new data
				this->OnImages(this->kinect_grabber_->GetLastDepth(), this->kinect_grabber_->GetLastColor());
			}
		});

}

void WdgtKinectSuperPoints::LoadRgbd(const std::string& fn)
{
	capture_filename_ = fn;

	StopProcessingThread();

	// load selected images
	std::string fn_color = fn + "_color.png";
	std::string fn_depth = fn + "_depth.pgm";
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

void WdgtKinectSuperPoints::StopProcessingThread()
{
	ui.horizontalSliderFrame->setEnabled(false);
	ui.horizontalSliderFrame->setMaximum(0);
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

void WdgtKinectSuperPoints::OnChangeFrame(int value)
{
	// if(mode_ == ReplayOniMode) {
	// 	this->reload = true;
	// }
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
	LoadRgbd(fn.toStdString());
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

void WdgtKinectSuperPoints::OnLive()
{
	ShowLive();
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
	if(ui.actionBatch_Save->isChecked()) {
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
	if(mode_ == LiveMode) {
		frame_counter_ ++;
	}
	has_new_frame_ = true;
	// save
	if(save_dasp_enabled_) {
		std::string fn1 = (boost::format(save_dasp_fn_+"%1$05d.tsv") % frame_counter_).str();
		SaveSuperpixels(dasp_processing_->clustering_, fn1, false);
		if(boost::num_vertices(dasp_processing_->Gnb) > 0) {
			std::string fn2 = (boost::format(save_dasp_fn_+"%1$05d_graph.txt") % frame_counter_).str();
			dasp::SaveGraph(dasp_processing_->Gnb, fn2);
		}
		else {
			if(frame_counter_ == 1) {
				std::cerr << "Graph not created -> not saving graph." << std::endl;
			}
		}
		{	// save rgb in png
			std::string fn_rgb = (boost::format(save_dasp_fn_+"%1$05d.png") % frame_counter_).str();
			slimage::Save(kinect_color, fn_rgb);
		}
	}
}

void WdgtKinectSuperPoints::OnUpdateImages()
{
	ui.action_Parameters->setChecked(gui_params_->isVisible());
	ui.action_Benchmarks->setChecked(gui_benchmark_->isVisible());
//	ui.actionAbout->setChecked(gui_benchmark_->isVisible());

	if(!has_new_frame_) {
		return;
	}

	if(mode_ == ReplayOniMode) {
		ui.horizontalSliderFrame->setValue(frame_counter_);
	}
	ui.labelFrameCounter->setText(QString("Frame: %1 / %2").arg(frame_counter_).arg(ui.horizontalSliderFrame->maximum()));
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

void WdgtKinectSuperPoints::onViewParameters()
{
	gui_params_->setVisible(ui.action_Parameters->isChecked());
}

void WdgtKinectSuperPoints::onViewBenchmark()
{
	gui_benchmark_->setVisible(ui.action_Benchmarks->isChecked());
}

void WdgtKinectSuperPoints::onViewAbout()
{
//	gui_about_->setVisible(ui.actionAbout->isChecked());
}

void WdgtKinectSuperPoints::closeEvent(QCloseEvent* event)
{
	gui_params_->close();
	gui_benchmark_->close();
//	gui_about_->close();
	event->accept();
}
