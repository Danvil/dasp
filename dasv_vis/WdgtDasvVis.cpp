#include "WdgtDasvVis.h"
#include <Slimage/IO.hpp>
#include <boost/bind.hpp>

void PrepareEngine(const boost::shared_ptr<Candy::Engine>& engine)
{
	// auto view_ = engine->getView();
	auto scene = engine->getScene();

	// boost::shared_ptr<Candy::DirectionalLight> light1(new Candy::DirectionalLight(Danvil::ctLinAlg::Vec3f(+1.0f, +1.0f, -1.0f)));
	// light1->setDiffuse(Danvil::Colorf(1.0f, 1.0f, 1.0f));
	// scene->addLight(light1);
	// boost::shared_ptr<Candy::DirectionalLight> light2(new Candy::DirectionalLight(Danvil::ctLinAlg::Vec3f(-1.0f, -1.0f, -1.0f)));
	// light2->setDiffuse(Danvil::Colorf(1.0f, 1.0f, 1.0f));
	// scene->addLight(light2);

	engine->setClearColor(Danvil::Color::Grey);

	scene->setShowCoordinateCross(true);

}

WdgtDasvVis::WdgtDasvVis(QWidget *parent)
    : QMainWindow(parent)
{
	ui.setupUi(this);

	std::cout << "Preparing OpenGL widgets ..." << std::endl;

	engine_main_ = ui.widgetCandy->getEngine();
	PrepareEngine(engine_main_);

	widget_candy_tab_ = new Candy::GLSystemQtWindow(0);
	ui.tabs->addTab(widget_candy_tab_, "3D");
	engine_tab_ = widget_candy_tab_->getEngine();
	PrepareEngine(engine_tab_);

	QObject::connect(&timer_tick_, SIGNAL(timeout()), this, SLOT(tick()));
	timer_tick_.setInterval(1);
	timer_tick_.start();

	dasv::DebugSetDisplayImageCallback(boost::bind(&WdgtDasvVis::showImage, this, _1, _2));
}

WdgtDasvVis::~WdgtDasvVis()
{
}

void WdgtDasvVis::setRgbdStream(const std::shared_ptr<RgbdStream>& rgbd_stream)
{
	rgbd_stream_ = rgbd_stream;
	dasv_ = std::make_shared<dasv::ContinuousSupervoxels>();
	dasv_->start();
}

void WdgtDasvVis::tick()
{
	if(rgbd_stream_->grab()) {
		Rgbd data = rgbd_stream_->get();
		showImage("color", data.color);
		// slimage::gui::Show("depth", data.depth, 500, 3000, 0);
		dasv_->step(data.color, data.depth);
	}
}

void WdgtDasvVis::showImage(const std::string& tag, const slimage::Image3ub& img)
{
	std::map<std::string, QWidget*> tabs_labels;
	// prepare tabs
	for(int i=0; i<ui.tabs->count(); i++) {
		std::string text = ui.tabs->tabText(i).toStdString();
		tabs_labels[text] = ui.tabs->widget(i);
	}
	// convert to Qt image
	QImage* qimg = slimage::qt::ConvertToQt(img);
	if(qimg == 0) {
		return;
	}
	// find label and create new if none exists
	QLabel* qlabel;
	if(tabs_labels.find(tag) != tabs_labels.end()) {
		qlabel = (QLabel*)tabs_labels[tag];
	} else {
		qlabel = new QLabel();
		tabs_labels[tag] = qlabel;
		ui.tabs->addTab(qlabel, QString::fromStdString(tag));
	}
	// display image in label
	qlabel->setPixmap(QPixmap::fromImage(*qimg));
	// cleanup
	delete qimg;
}
