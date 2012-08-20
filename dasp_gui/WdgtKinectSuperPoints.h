#ifndef WDGTKINECTSUPERPOINTS_H
#define WDGTKINECTSUPERPOINTS_H

#include "ui_WdgtKinectSuperPoints.h"
#include "WdgtSuperpixelParameters.h"
#include "WdgtBenchmark.h"
#if defined DASP_HAS_OPENNI
#	include "KinectGrabber.h"
#endif
#include "DaspProcessing.h"
#if defined DASP_HAS_SIMPLEENGINE
#	include <Danvil/SimpleEngine/System/GLSystemQtWindow.h>
#	include <Danvil/SimpleEngine.h>
#endif
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
	void OnImages(const slimage::Image1ui16& kinect_depth, const slimage::Image3ub& kinect_color);

	void ComputeBlueNoiseImpl();

public Q_SLOTS:
	void OnUpdateImages();
	void OnLoadOne();
#if defined DASP_HAS_OPENNI
	void OnCaptureOne();
	void OnLoadOni();
	void OnLive();
#endif
	void OnSaveDebugImages();

private:
#if defined DASP_HAS_SIMPLEENGINE
 	PTR(Danvil::SimpleEngine::View) view_;
	PTR(Danvil::SimpleEngine::Scene) scene_;
	PTR(Danvil::SimpleEngine::Engine) engine_;
	Danvil::SimpleEngine::GLSystemQtWindow* gl_wdgt_;
#endif

	QTimer timer_;

	boost::shared_ptr<WdgtSuperpixelParameters> gui_params_;
	boost::shared_ptr<WdgtBenchmark> gui_benchmark_;

#if defined DASP_HAS_OPENNI
	boost::shared_ptr<Romeo::Kinect::KinectGrabber> kinect_grabber_;
#endif

	boost::shared_ptr<DaspProcessing> dasp_processing_;

	QMutex images_mutex_;

	boost::thread kinect_thread_;

	bool capture_next_;
	std::string capture_filename_;
	bool interrupt_loaded_thread_;
	bool save_debug_next_;
	bool reload;

private:
    Ui::WdgtKinectSuperPointsClass ui;
};

#endif // WDGTKINECTSUPERPOINTS_H
