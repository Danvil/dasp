#ifndef WDGTKINECTCONTROL_H
#define WDGTKINECTCONTROL_H

#include <QtGui/QWidget>
#include "ui_WdgtSuperpixelParameters.h"
#include <SuperPoints/Superpixels.hpp>

class WdgtSuperpixelParameters : public QWidget
{
    Q_OBJECT

public:
    WdgtSuperpixelParameters(const boost::shared_ptr<dasp::Parameters>& dasp_params, QWidget *parent = 0);
    ~WdgtSuperpixelParameters();

public Q_SLOTS:
	void ChangeSuperpixelCount(int val);
	void ChangeSuperpixelIterations(int val);
	void ChangeSuperpixelWeightSpatial(double val);
	void ChangeSuperpixelWeightDepth(double val);
	void ChangeSuperpixelWeightNormal(double val);
	void ChangeSuperpixelCoverage(double val);
	void ChangeSuperpixelDensity(int val);

private:
    boost::shared_ptr<dasp::Parameters> dasp_params_;

private:
    Ui::WdgtSuperpixelParametersClass ui;
};

#endif
