#ifndef WDGTKINECTCONTROL_H
#define WDGTKINECTCONTROL_H

#include <QtGui/QWidget>
#include "ui_WdgtSuperpixelParameters.h"
#include <SuperPoints/Superpixels.hpp>
#include <functional>

class WdgtSuperpixelParameters : public QWidget
{
    Q_OBJECT

public:
    WdgtSuperpixelParameters(const boost::shared_ptr<dasp::Parameters>& dasp_params, QWidget *parent = 0);
    ~WdgtSuperpixelParameters();

    std::function<void()> on_train_;
    std::function<void(float)> on_change_cm_sigma_scale_;

public Q_SLOTS:
	void OnSuperSeedType(const QString& txt);
	void ChangeSuperpixelCountOrDensity(int val);
	void ChangeSuperpixelIterations(int val);
	void ChangeSuperpixelWeightColor(double val);
	void ChangeSuperpixelWeightSpatial(double val);
	void ChangeSuperpixelWeightDepth(double val);
	void ChangeSuperpixelWeightNormal(double val);
	void ChangeSuperpixelCoverage(double val);
	void ChangeColorModelSigmaScale(double val);
	void OnColorModelTrain();

private:
    boost::shared_ptr<dasp::Parameters> dasp_params_;

private:
    Ui::WdgtSuperpixelParametersClass ui;
};

#endif
