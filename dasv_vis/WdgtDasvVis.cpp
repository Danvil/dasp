#include "WdgtDasvVis.h"
#include <graphseg/Rendering.hpp>
#include <Slimage/IO.hpp>
#include <QtGui/QMdiSubWindow>
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

	widget_candy_global_ = new Candy::GLSystemQtWindow(0);
	auto w1 = ui.mdiArea->addSubWindow(widget_candy_global_);
	w1->setWindowTitle("3D global");
	engine_global_ = widget_candy_global_->getEngine();
	PrepareEngine(engine_global_);

	widget_candy_frame_ = new Candy::GLSystemQtWindow(0);
	auto w2 = ui.mdiArea->addSubWindow(widget_candy_frame_);
	w2->setWindowTitle("3D frame");
	engine_frame_ = widget_candy_frame_->getEngine();
	PrepareEngine(engine_frame_);

	{
		boost::shared_ptr<Candy::IRenderable> renderling(new Candy::ObjectRenderling(
			[this]() {
				this->renderGraphGlobal();
			}
		));
		engine_global_->getScene()->addItem(renderling);
	}

	{
		boost::shared_ptr<Candy::IRenderable> renderling(new Candy::ObjectRenderling(
			[this]() {
				this->renderGraphFrame();
			}
		));
		engine_frame_->getScene()->addItem(renderling);
	}

	QObject::connect(&timer_tick_, SIGNAL(timeout()), this, SLOT(tick()));
	timer_tick_.setInterval(50);
	timer_tick_.start();

	dasv::DebugSetDisplayImageCallback(boost::bind(&WdgtDasvVis::showImageThreadsafe, this, _1, _2));

}

WdgtDasvVis::~WdgtDasvVis()
{
}

void WdgtDasvVis::setRgbdStream(const std::shared_ptr<RgbdStream>& rgbd_stream)
{
	rgbd_stream_ = rgbd_stream;

	worker_interupt_ = true;
	if(worker_.joinable())
		worker_.join();

	worker_ = std::thread(
		[this]() {
			worker_interupt_ = false;
			this->dasv_ = std::make_shared<dasv::ContinuousSupervoxels>();
			this->dasv_->start();
			while(!worker_interupt_ && this->rgbd_stream_->grab()) {
				this->dasv_->step(this->rgbd_stream_->get());
				{
					std::lock_guard<std::mutex> lock(dasv_graph_mutex_);
					dasv_graph_ = dasv_->getGraph();
				}
			}
		});

}

void WdgtDasvVis::tick()
{
	std::lock_guard<std::mutex> lock(show_images_cache_mutex_);
	if(show_images_cache_.size() == 0)
		return;
	for(auto p : show_images_cache_) {
		showImage(p.first, p.second);
	}
	std::cout << "Updated " << show_images_cache_.size() << " images" << std::endl;
	show_images_cache_.clear();
}

void WdgtDasvVis::showImageThreadsafe(const std::string& tag, const slimage::Image3ub& img)
{
	std::lock_guard<std::mutex> lock(show_images_cache_mutex_);
	show_images_cache_[tag] = img;
}

void WdgtDasvVis::showImage(const std::string& tag, const slimage::Image3ub& img)
{
	// convert to Qt image
	QImage* qimg = slimage::qt::ConvertToQt(img);
	if(qimg == 0) {
		return;
	}
	// prepare subwindow list
	std::map<std::string, QMdiSubWindow*> subwindows_by_tag;
	for(QMdiSubWindow* w : ui.mdiArea->subWindowList()) {
		subwindows_by_tag[w->windowTitle().toStdString()] = w;
	}
	// find label
	auto it = subwindows_by_tag.find(tag);
	// create new if none exists
	QMdiSubWindow* sw;
	QLabel* qlabel;
	if(it != subwindows_by_tag.end()) {
		sw = it->second;
		qlabel = (QLabel*)sw->widget();
	} else {
		qlabel = new QLabel();
		subwindows_by_tag[tag] = sw;
		sw = ui.mdiArea->addSubWindow(qlabel);
		sw->setWindowTitle(QString::fromStdString(tag));
		sw->show();
	}
	// display image in label
	qlabel->setPixmap(QPixmap::fromImage(*qimg));
	qlabel->adjustSize();
	sw->adjustSize();
	// cleanup
	delete qimg;
}

void WdgtDasvVis::renderGraphGlobal()
{
	if(!dasv_) {
		return;
	}
	dasv::ClusterGraph graph;
	{
		std::lock_guard<std::mutex> lock(dasv_graph_mutex_);
		graph = dasv_graph_;
	}
	graphseg::RenderEdges3DVcol(graph,
		[&graph](const dasv::ClusterGraph::vertex_descriptor& vid) {
			// returns vertex coordinate
			const dasv::Cluster& c = graph[vid];
			return Eigen::Vector3f{ 0.01f*c.pixel.x(), 0.01f*c.pixel.y(), 0.03f*static_cast<float>(c.time) };
			// return c.position;
		},
		[&graph](const dasv::ClusterGraph::vertex_descriptor& vid) {
			// returns vertex color
			return graph[vid].color;
		}
	);

}

void WdgtDasvVis::renderGraphFrame()
{
	if(!dasv_) {
		return;
	}
	dasv::ClusterGraph graph;
	{
		std::lock_guard<std::mutex> lock(dasv_graph_mutex_);
		graph = dasv_graph_;
	}
	graphseg::RenderEdges3DVcol(graph,
		[&graph](const dasv::ClusterGraph::vertex_descriptor& vid) {
			// returns vertex coordinate
			return graph[vid].position;
		},
		[&graph](const dasv::ClusterGraph::vertex_descriptor& vid) {
			// returns vertex color
			return graph[vid].color;
		}
	);

}
