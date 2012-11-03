#ifndef WDGTKINECTCONTROL_H
#define WDGTKINECTCONTROL_H

#include <QtGui/QWidget>
#include "ui_WdgtSuperpixelParameters.h"
#include "DaspProcessing.h"
#include <functional>

class WdgtSuperpixelParameters : public QWidget
{
    Q_OBJECT

public:
	WdgtSuperpixelParameters(const boost::shared_ptr<DaspProcessing>& dasp_processing, QWidget *parent = 0);
	~WdgtSuperpixelParameters();

	bool* reload;

	void SetActualCount(unsigned int count);
	void SetActualRadius(float radius);

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
	void ChangePlotGraphSpatialCut(int state);
	void ChangePlotGraph(int state);
	void ChangePlotGraphWeights(int state);
	void ChangePlotDensity(int state);
	void ChangePlotSegments(int state);
	void ChangeClipEnable(int state);
	void ChangeClipXMin(double val);
	void ChangeClipYMin(double val);
	void ChangeClipZMin(double val);
	void ChangeClipXMax(double val);
	void ChangeClipYMax(double val);
	void ChangeClipZMax(double val);

private:
	boost::shared_ptr<DaspProcessing> dasp_processing_;

private:
    Ui::WdgtSuperpixelParametersClass ui;
};

#endif
