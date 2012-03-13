#ifndef WDGTKINECTSUPERPOINTS_H
#define WDGTKINECTSUPERPOINTS_H

#include "ui_WdgtKinectSuperPoints.h"
#include "WdgtSuperpixelParameters.h"
#include <SuperPoints/DaspTracker.h>
#include <Romeo/Kinect/KinectGrabber.h>
#include <Danvil/SimpleEngine/System/GLSystemQtWindow.h>
#include <Danvil/SimpleEngine.h>
#include <QtGui/QMainWindow>
#include <QtCore/QTimer>
#include <QtCore/QMutex>
#include <boost/thread.hpp>

class WdgtKinectSuperPoints : public QMainWindow
{
    Q_OBJECT

public:
    WdgtKinectSuperPoints(QWidget *parent = 0);
    ~WdgtKinectSuperPoints();

private:
	void OnImages(Danvil::Images::Image1ui16Ptr kinect_depth, Danvil::Images::Image3ubPtr kinect_color);

	void ComputeBlueNoiseImpl();

public Q_SLOTS:
	void OnUpdateImages();
	void OnCaptureOne();
	void OnLoadOne();
	void OnLoadOni();
	void OnLive();
	void OnSaveDebugImages();

private:
 	PTR(Danvil::SimpleEngine::View) view_;
	PTR(Danvil::SimpleEngine::Scene) scene_;
	PTR(Danvil::SimpleEngine::Engine) engine_;
	Danvil::SimpleEngine::GLSystemQtWindow* gl_wdgt_;

	QTimer timer_;

	boost::shared_ptr<WdgtSuperpixelParameters> gui_params_;
	boost::shared_ptr<Romeo::Kinect::KinectGrabber> kinect_grabber_;

	boost::shared_ptr<dasp::DaspTracker> dasp_tracker_;

	QMutex images_mutex_;

	boost::thread kinect_thread_;

	bool capture_next_;
	std::string capture_filename_;
	bool interrupt_loaded_thread_;
	bool save_debug_next_;

private:
    Ui::WdgtKinectSuperPointsClass ui;
};

#endif // WDGTKINECTSUPERPOINTS_H
