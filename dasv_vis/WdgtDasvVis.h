#ifndef DASVVIS_WDGTDASVVIS_H
#define DASVVIS_WDGTDASVVIS_H

#include "ui_WdgtDasvVis.h"
#include <dasv.hpp>
#include <rgbd.hpp>
#include <Candy/System/GLSystemQtWindow.h>
#include <Candy.h>
#include <Slimage/Slimage.hpp>
#include <QtGui/QMainWindow>
#include <QtCore/QTimer>
#include <thread>
#include <mutex>

class WdgtDasvVis : public QMainWindow
{
    Q_OBJECT

public:
	WdgtDasvVis(QWidget *parent = 0);
	~WdgtDasvVis();

public:
	void setRgbdStream(const std::shared_ptr<RgbdStream>& rgbd_stream);

public Q_SLOTS:
	void tick();

private:
	void showImageThreadsafe(const std::string& tag, const slimage::Image3ub& img);
	void showImage(const std::string& tag, const slimage::Image3ub& img);
	void renderGraphGlobal();
	void renderGraphFrame();

private:
	Candy::GLSystemQtWindow* widget_candy_global_;
	boost::shared_ptr<Candy::Engine> engine_global_;
	Candy::GLSystemQtWindow* widget_candy_frame_;
	boost::shared_ptr<Candy::Engine> engine_frame_;

	std::shared_ptr<RgbdStream> rgbd_stream_;

	std::shared_ptr<dasv::ContinuousSupervoxels> dasv_;

	QTimer timer_tick_;

	std::thread worker_;
	bool worker_interupt_;
	std::map<std::string, slimage::Image3ub> show_images_cache_;
	std::mutex show_images_cache_mutex_;

	dasv::ClusterGraph dasv_graph_;
	std::mutex dasv_graph_mutex_;

private:
    Ui::WdgtDasvVisClass ui;
};

#endif
