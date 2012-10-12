#ifndef WDGTKINECTSUPERPOINTS_H
#define WDGTKINECTSUPERPOINTS_H

#include "ui_WdgtKinectSuperPoints.h"
#include "WdgtSuperpixelParameters.h"
#include "WdgtBenchmark.h"
#if defined DASP_HAS_OPENNI
#	include "KinectGrabber.h"
#endif
#include "DaspProcessing.h"
#if defined DASP_HAS_CANDY
#	include <Candy/System/GLSystemQtWindow.h>
#	include <Candy.h>
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

	void StopProcessingThread();

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
#if defined DASP_HAS_CANDY
 	PTR(Candy::View) view_;
	PTR(Candy::Scene) scene_;
	PTR(Candy::Engine) engine_;
	Candy::GLSystemQtWindow* gl_wdgt_;
#endif

	QTimer timer_;

	boost::shared_ptr<WdgtSuperpixelParameters> gui_params_;
	boost::shared_ptr<WdgtBenchmark> gui_benchmark_;

#if defined DASP_HAS_OPENNI
	boost::shared_ptr<Romeo::Kinect::KinectGrabber> kinect_grabber_;
#endif

	boost::shared_ptr<DaspProcessing> dasp_processing_;

	boost::thread processing_thread_;

	bool capture_next_;
	std::string capture_filename_;
	bool interrupt_loaded_thread_;
	bool reload;

	unsigned int frame_counter_;
	bool has_new_frame_;

private:
    Ui::WdgtKinectSuperPointsClass ui;
};

#endif // WDGTKINECTSUPERPOINTS_H
