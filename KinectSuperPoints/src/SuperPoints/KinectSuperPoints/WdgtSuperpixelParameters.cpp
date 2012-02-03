#include "WdgtSuperpixelParameters.h"

WdgtSuperpixelParameters::WdgtSuperpixelParameters(const boost::shared_ptr<dasp::Parameters>& dasp_params, QWidget *parent)
    : QWidget(parent)
{
	dasp_params_ = dasp_params;

	ui.setupUi(this);

	QObject::connect(ui.spinBoxCount, SIGNAL(valueChanged(int)), this, SLOT(ChangeSuperpixelCount(int)));
	QObject::connect(ui.spinBoxIterations, SIGNAL(valueChanged(int)), this, SLOT(ChangeSuperpixelIterations(int)));
	QObject::connect(ui.doubleSpinBoxWeightSpatial, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightSpatial(double)));
	QObject::connect(ui.doubleSpinBoxWeightDepth, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightDepth(double)));
	QObject::connect(ui.doubleSpinBoxWeightNormal, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightNormal(double)));
	QObject::connect(ui.doubleSpinBoxCoverage, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelCoverage(double)));
	QObject::connect(ui.spinBoxDensity, SIGNAL(valueChanged(int)), this, SLOT(ChangeSuperpixelDensity(int)));

	dasp_params_->cluster_count = ui.spinBoxCount->value();
	dasp_params_->iterations = ui.spinBoxIterations->value();
	dasp_params_->weight_spatial = ui.doubleSpinBoxWeightSpatial->value();
	dasp_params_->weight_depth = ui.doubleSpinBoxWeightDepth->value();
	dasp_params_->weight_normal = ui.doubleSpinBoxWeightNormal->value();
	dasp_params_->coverage = ui.doubleSpinBoxCoverage->value();
	dasp_params_->roh = ui.spinBoxDensity->value();

}

WdgtSuperpixelParameters::~WdgtSuperpixelParameters()
{

}

void WdgtSuperpixelParameters::ChangeSuperpixelCount(int val)
{
	dasp_params_->cluster_count = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelIterations(int val)
{
	dasp_params_->iterations = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightSpatial(double val)
{
	dasp_params_->weight_spatial = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightDepth(double val)
{
	dasp_params_->weight_depth = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightNormal(double val)
{
	dasp_params_->weight_normal = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelCoverage(double val)
{
	dasp_params_->coverage = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelDensity(int val)
{
	dasp_params_->roh = val;
}
