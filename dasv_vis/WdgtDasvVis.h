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
	void showImage(const std::string& tag, const slimage::Image3ub& img);

private:
 	boost::shared_ptr<Candy::Engine> engine_main_;
	boost::shared_ptr<Candy::Engine> engine_tab_;
	Candy::GLSystemQtWindow* widget_candy_tab_;

	std::shared_ptr<RgbdStream> rgbd_stream_;

	std::shared_ptr<dasv::ContinuousSupervoxels> dasv_;

	QTimer timer_tick_;

private:
    Ui::WdgtDasvVisClass ui;
};

#endif
