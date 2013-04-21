#include "WdgtKinectSuperPoints.h"
#include <dasp/IO.hpp>
#include <graphseg/IO.hpp>
#include <Slimage/IO.hpp>
#include <Slimage/Parallel.h>
#include <boost/bind.hpp>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <boost/format.hpp>

WdgtKinectSuperPoints::WdgtKinectSuperPoints(bool no3d, QWidget *parent)
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
	QObject::connect(ui.actionSave_Active_Image, SIGNAL(triggered()), this, SLOT(OnSaveActiveImage()));
	QObject::connect(ui.actionBatch_Save, SIGNAL(triggered()), this, SLOT(OnSaveDasp()));

	QObject::connect(ui.action_Settings, SIGNAL(triggered()), this, SLOT(onViewSettings()));
	QObject::connect(ui.action_DaspParameters, SIGNAL(triggered()), this, SLOT(onViewDaspParameters()));
	QObject::connect(ui.action_Benchmarks, SIGNAL(triggered()), this, SLOT(onViewBenchmark()));
	QObject::connect(ui.actionAbout, SIGNAL(triggered()), this, SLOT(onViewAbout()));

	dasp_processing_.reset(new DaspProcessing());

	gui_dasp_params_.reset(new WdgtDaspParameters(dasp_processing_->dasp_params));
	gui_dasp_params_->setAttribute(Qt::WA_DeleteOnClose, false);
	gui_dasp_params_->reload = &reload;
	gui_dasp_params_->show();

	gui_settings_.reset(new WdgtSettings(dasp_processing_));
	gui_settings_->setAttribute(Qt::WA_DeleteOnClose, false);
	gui_settings_->reload = &reload;
	gui_settings_->show();

	gui_benchmark_.reset(new WdgtBenchmark());
	gui_benchmark_->setAttribute(Qt::WA_DeleteOnClose, false);
	Danvil::Benchmark::Instance().setOnUpdate(boost::bind(&WdgtBenchmark::update, gui_benchmark_, _1, _2));

	gui_about_.reset(new WdgtAbout());
	gui_about_->setAttribute(Qt::WA_DeleteOnClose, false);

#if defined DASP_HAS_CANDY
	if(!no3d) {
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
	}
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
	std::string config = "/home/david/git/OpenNI/Platform/Linux/Redist/OpenNI-Bin-Dev-Linux-x64-v1.5.4.0/Samples/Config/SamplesConfig.xml";
	std::cout << "Opening kinect config file: " << config << std::endl;
	rgbd_stream_ = FactorKinectLive(config);
	mode_ = LiveMode;

	// processing thread polls kinect
	frame_counter_ = 0;
	processing_thread_ = boost::thread(
		[this]() {
		this->interrupt_loaded_thread_ = false;
		while(!this->interrupt_loaded_thread_ && this->rgbd_stream_->grab()) {
				this->OnImages(this->rgbd_stream_->get());
			}
		});
}

void WdgtKinectSuperPoints::LoadOni(const std::string& fn)
{
	StopProcessingThread();

	// open Kinect OpenNI ONI file
	std::cout << "Opening oni file: " << fn << std::endl;
	auto random_access_rgbd_stream = FactorOni(fn);
	rgbd_stream_ = random_access_rgbd_stream;

	mode_ = ReplayOniMode;
	remembered_capture_fn_ = fn;
	
	ui.horizontalSliderFrame->setEnabled(true);
	std::cout << "Number of frames: " << random_access_rgbd_stream->numFrames() << std::endl;
	ui.horizontalSliderFrame->setMaximum(random_access_rgbd_stream->numFrames());

	// processing thread polls kinect
	frame_counter_ = 0;
	processing_thread_ = boost::thread(
		[this, random_access_rgbd_stream]() {
			this->interrupt_loaded_thread_ = false;
			while(!this->interrupt_loaded_thread_) {
				unsigned int requested_frame = frame_counter_ + (this->ui.checkBoxPlay->isChecked() ? 1 : 0);
				while(!this->reload && requested_frame == frame_counter_ && !this->interrupt_loaded_thread_) {
					usleep(1000);
					requested_frame = ui.horizontalSliderFrame->value() + (this->ui.checkBoxPlay->isChecked() ? 1 : 0);
				}
				this->reload = false;
				// seek to requested frame (playing will start from there)
				random_access_rgbd_stream->seek(requested_frame);
				frame_counter_ = random_access_rgbd_stream->tell();
				// grab frame from oni
				if(!this->rgbd_stream_->grab()) {
					break;
				}
				// update images with new data
				this->OnImages(this->rgbd_stream_->get());
			}
		});

}

void WdgtKinectSuperPoints::LoadRgbd(const std::string& fn)
{
	capture_filename_ = fn;

	StopProcessingThread();

	// load selected images
	std::cout << "Reading '" << fn << "_color.png' and '" << fn << "_depth.pgm'..." << std::flush;
	rgbd_stream_ = FactorStatic(fn);

	mode_ = SingleFileMode;

	// processing thread repeats loaded image (only updated if requested)
	processing_thread_ = boost::thread(
		[this]() {
			this->interrupt_loaded_thread_ = false;
			while(!this->interrupt_loaded_thread_) {
				this->OnImages(rgbd_stream_->get());
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
	interrupt_loaded_thread_ = true;
	processing_thread_.join();
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
	QString fn = QFileDialog::getSaveFileName(this, "Open Color/Depth files");
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

void WdgtKinectSuperPoints::OnSaveActiveImage()
{
	// ask for filename
	QString fn = QFileDialog::getSaveFileName(this,
		"Save Active Image", "/home/david");
	// save active generated dasp images
	((QLabel*)ui.tabs->widget(ui.tabs->currentIndex()))->pixmap()->save(fn);
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

void WdgtKinectSuperPoints::OnImages(const Rgbd& rgbd)
{
	// save RGBD image if requested
	if(capture_next_) {
		slimage::Save(rgbd.color, capture_filename_ + "_color.png");
		slimage::Save(rgbd.depth, capture_filename_ + "_depth.pgm");
		capture_next_ = false;
	}
	// process image
	dasp_processing_->step(rgbd.depth, rgbd.color);
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
			graphseg::WriteEdges(fn2,
				dasp_processing_->Gnb_weighted,
				boost::get(boost::edge_bundle, dasp_processing_->Gnb_weighted));
		}
		else {
			if(frame_counter_ == 1) {
				std::cerr << "Graph not created -> not saving graph." << std::endl;
			}
		}
		{	// save rgb in png
			std::string fn_rgb = (boost::format(save_dasp_fn_+"%1$05d.png") % frame_counter_).str();
			slimage::Save(rgbd.color, fn_rgb);
		}
	}
}

void WdgtKinectSuperPoints::OnUpdateImages()
{
	ui.action_DaspParameters->setChecked(gui_dasp_params_->isVisible());
	ui.action_Settings->setChecked(gui_settings_->isVisible());
	ui.action_Benchmarks->setChecked(gui_benchmark_->isVisible());
	ui.actionAbout->setChecked(gui_benchmark_->isVisible());

	if(!has_new_frame_) {
		return;
	}

	if(dasp_processing_->dasp_params->count == 0) {
		unsigned int actual_count = dasp_processing_->clustering_.cluster.size();
		gui_dasp_params_->SetActualCount(actual_count);
	}
	else {
		double actual_radius = dasp_processing_->clustering_.opt.base_radius;
		gui_dasp_params_->SetActualRadius(actual_radius);
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

void WdgtKinectSuperPoints::onViewDaspParameters()
{
	gui_dasp_params_->setVisible(ui.action_DaspParameters->isChecked());
}

void WdgtKinectSuperPoints::onViewSettings()
{
	gui_settings_->setVisible(ui.action_Settings->isChecked());
}

void WdgtKinectSuperPoints::onViewBenchmark()
{
	gui_benchmark_->setVisible(ui.action_Benchmarks->isChecked());
}

void WdgtKinectSuperPoints::onViewAbout()
{
	gui_about_->setVisible(ui.actionAbout->isChecked());
}

void WdgtKinectSuperPoints::closeEvent(QCloseEvent* event)
{
	gui_dasp_params_->close();
	gui_settings_->close();
	gui_benchmark_->close();
	gui_about_->close();
	event->accept();
}
