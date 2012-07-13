#ifndef WDGTKINECTCONTROL_H
#define WDGTKINECTCONTROL_H

#include <QtGui/QWidget>
#include "ui_WdgtSuperpixelParameters.h"
#include <dasp/DaspTracker.h>
#include <functional>

class WdgtSuperpixelParameters : public QWidget
{
    Q_OBJECT

public:
	WdgtSuperpixelParameters(const boost::shared_ptr<dasp::DaspTracker>& dasp_tracker, QWidget *parent = 0);
	~WdgtSuperpixelParameters();

	bool* reload;

public Q_SLOTS:
	void ChangeDaspSmoothDepth(int state);
	void ChangeDaspRepairDepth(int state);
	void OnSuperSeedType(int selection);
	void ChangeSuperUseGradientDensity(int state);
	void ChangeSuperpixelSkipBad(int state);
	void ChangeSuperpixelRadius(double val);
	void ChangeSuperpixelCount(int val);
	void ChangeSuperpixelIterations(int val);
	void ChangeSuperpixelWeightSpatial(double val);
	void ChangeSuperpixelWeightColor(double val);
	void OnDaspColorSpace(int selection);
	void ChangeSuperpixelWeightNormal(double val);
	void ChangeSuperpixelWeightDepth(double val);
	void ChangeSuperConquerEnclaves(int val);
	void ChangeSuperpixelCoverage(double val);
	void ChangeDaspSegmentThreshold(double val);
	void ChangePlotPoints(int state);
	void ChangePlotPointsColor(int selection);
	void ChangePlotClusters(int state);
	void ChangePlotClusterMode(int selection);
	void ChangePlotClusterColor(int selection);
	void ChangePlotBorders(int state);
	void ChangePlotGraph(int state);
	void ChangePlotDensity(int state);
	void ChangePlotSegments(int state);

private:
	boost::shared_ptr<dasp::DaspTracker> dasp_tracker_;

private:
    Ui::WdgtSuperpixelParametersClass ui;
};

#endif
