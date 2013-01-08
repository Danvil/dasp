#ifndef WDGTDASPPARAMETERS_H
#define WDGTDASPPARAMETERS_H

#include <QtGui/QWidget>
#include "ui_WdgtDaspParameters.h"
#include "dasp/Parameters.hpp"
#include <boost/shared_ptr.hpp>

class WdgtDaspParameters : public QWidget
{
    Q_OBJECT

public:
	WdgtDaspParameters(const boost::shared_ptr<dasp::Parameters>& dasp_opt, QWidget *parent = 0);
	~WdgtDaspParameters();

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
	void ChangeClipEnable(int state);
	void ChangeClipXMin(double val);
	void ChangeClipYMin(double val);
	void ChangeClipZMin(double val);
	void ChangeClipXMax(double val);
	void ChangeClipYMax(double val);
	void ChangeClipZMax(double val);

private:
	boost::shared_ptr<dasp::Parameters> dasp_opt_;

private:
    Ui::WdgtDaspParametersClass ui;
};

#endif
