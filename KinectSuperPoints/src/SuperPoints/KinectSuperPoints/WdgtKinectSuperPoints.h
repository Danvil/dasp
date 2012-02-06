#ifndef WDGTKINECTSUPERPOINTS_H
#define WDGTKINECTSUPERPOINTS_H

#include "ui_WdgtKinectSuperPoints.h"
#include "WdgtSuperpixelParameters.h"
#include <Romeo/Kinect/KinectGrabber.h>
#include <Slimage/Slimage.hpp>
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

private:
	QTimer timer_;

	boost::shared_ptr<WdgtSuperpixelParameters> gui_params_;
    boost::shared_ptr<dasp::Parameters> dasp_params;
	boost::shared_ptr<Romeo::Kinect::KinectGrabber> kinect_grabber_;

	std::map<std::string, slimage::ImagePtr> images_;

	QMutex images_mutex_;

	boost::thread kinect_thread_;
	bool running_;

private:
    Ui::WdgtKinectSuperPointsClass ui;
};

#endif // WDGTKINECTSUPERPOINTS_H
