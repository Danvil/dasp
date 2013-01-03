#include "WdgtDasvVis.h"

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
}

WdgtDasvVis::~WdgtDasvVis()
{
}
