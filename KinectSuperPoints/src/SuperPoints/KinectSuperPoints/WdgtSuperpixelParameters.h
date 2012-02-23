#ifndef WDGTKINECTCONTROL_H
#define WDGTKINECTCONTROL_H

#include <QtGui/QWidget>
#include "ui_WdgtSuperpixelParameters.h"
#include <SuperPoints/DaspTracker.h>
#include <functional>

class WdgtSuperpixelParameters : public QWidget
{
    Q_OBJECT

public:
	WdgtSuperpixelParameters(const boost::shared_ptr<dasp::DaspTracker>& dasp_tracker, QWidget *parent = 0);
	~WdgtSuperpixelParameters();

public Q_SLOTS:
	void OnSuperSeedType(int selection);
	void ChangeSuperUseGradientDensity(int state);
	void ChangeSuperpixelRadius(double val);
	void ChangeSuperpixelIterations(int val);
	void ChangeSuperpixelWeightColor(double val);
	void ChangeSuperpixelWeightSpatial(double val);
	void ChangeSuperpixelWeightDepth(double val);
	void ChangeSuperpixelWeightNormal(double val);
	void ChangeSuperpixelCoverage(double val);
	void ChangeColorModelSigmaScale(double val);
	void OnColorModelTrain();
	void ChangePlotPoints(int state);
	void ChangePlotPointsColor(int selection);
	void ChangePlotClusters(int state);
	void ChangePlotClusterMode(int selection);
	void ChangePlotClusterColor(int selection);
	void ChangePlotBorders(int state);
	void ChangePlotGraph(int state);

private:
	boost::shared_ptr<dasp::DaspTracker> dasp_tracker_;

private:
    Ui::WdgtSuperpixelParametersClass ui;
};

#endif
