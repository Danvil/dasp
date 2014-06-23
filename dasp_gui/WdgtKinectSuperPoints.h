#ifndef WDGTKINECTSUPERPOINTS_H
#define WDGTKINECTSUPERPOINTS_H

#include "ui_WdgtKinectSuperPoints.h"
#include "WdgtSettings.h"
#include "WdgtBenchmark.h"
#include "WdgtAbout.h"
#include "DaspProcessing.h"
#include <common/WdgtDaspParameters.h>
#include <rgbd.hpp>
#if defined DASP_HAS_CANDY
#	include <Candy/System/GLSystemQtWindow.h>
#	include <Candy.h>
#endif
#include <QtGui/QMainWindow>
#include <QtGui/QCloseEvent>
#include <QtCore/QTimer>
#include <QtCore/QMutex>
#include <boost/thread.hpp>
#include <memory>

class WdgtKinectSuperPoints : public QMainWindow
{
    Q_OBJECT

public:
	WdgtKinectSuperPoints(bool no3d, QWidget *parent = 0);
	~WdgtKinectSuperPoints();

	void ShowLive();

	void LoadOni(const std::string& fn);

	void LoadRgbd(const std::string& fn);

private:
	void OnImages(const Rgbd& rgbd);

	void StopProcessingThread();

public Q_SLOTS:
	void OnChangeFrame(int);
	void OnUpdateImages();
	void OnLoadOne();
#if defined DASP_HAS_OPENNI
	void OnCaptureOne();
	void OnLoadOni();
	void OnLive();
#endif
	void OnSaveDebugImages();
	void OnSaveActiveImage();
	void OnSaveSuperpixels();
	void OnSaveDasp();
	void onViewDaspParameters();
	void onViewSettings();
	void onViewBenchmark();
	void onViewAbout();

protected:
	void closeEvent(QCloseEvent* event);

private:
#if defined DASP_HAS_CANDY
 	PTR(Candy::View) view_;
	PTR(Candy::Scene) scene_;
	PTR(Candy::Engine) engine_;
	Candy::GLSystemQtWindow* gl_wdgt_;
#endif

	QTimer timer_;

	boost::shared_ptr<WdgtDaspParameters> gui_dasp_params_;
	boost::shared_ptr<WdgtSettings> gui_settings_;
	boost::shared_ptr<WdgtBenchmark> gui_benchmark_;
	boost::shared_ptr<WdgtAbout> gui_about_;

	std::shared_ptr<RgbdStream> rgbd_stream_;

	boost::shared_ptr<DaspProcessing> dasp_processing_;

	boost::thread processing_thread_;

	enum CaptureMode {
		IdleMode, SingleFileMode, ReplayOniMode, LiveMode
	};
	CaptureMode mode_;

	std::string remembered_capture_fn_;

	bool capture_next_;
	std::string capture_filename_;
	bool interrupt_loaded_thread_;
	bool reload;

	unsigned int frame_counter_;
	bool has_new_frame_;
	bool save_dasp_enabled_;
	std::string save_dasp_fn_;

private:
    Ui::WdgtKinectSuperPointsClass ui;
};

#endif // WDGTKINECTSUPERPOINTS_H
