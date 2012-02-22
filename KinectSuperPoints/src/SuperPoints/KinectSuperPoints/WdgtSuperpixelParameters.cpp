#include "WdgtSuperpixelParameters.h"

WdgtSuperpixelParameters::WdgtSuperpixelParameters(const boost::shared_ptr<dasp::Parameters>& dasp_params, QWidget *parent)
    : QWidget(parent)
{
	dasp_params_ = dasp_params;

	ui.setupUi(this);

	QObject::connect(ui.comboBoxSeedType, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(OnSuperSeedType(const QString&)));
	QObject::connect(ui.doubleSpinBoxRadius, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelRadius(double)));
	QObject::connect(ui.spinBoxIterations, SIGNAL(valueChanged(int)), this, SLOT(ChangeSuperpixelIterations(int)));
	QObject::connect(ui.doubleSpinBoxWeightColor, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightColor(double)));
	QObject::connect(ui.doubleSpinBoxWeightSpatial, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightSpatial(double)));
	QObject::connect(ui.doubleSpinBoxWeightDepth, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightDepth(double)));
	QObject::connect(ui.doubleSpinBoxWeightNormal, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightNormal(double)));
	QObject::connect(ui.doubleSpinBoxCoverage, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelCoverage(double)));
	QObject::connect(ui.pushButtonColorModelTrain, SIGNAL(clicked()), this, SLOT(OnColorModelTrain()));
	QObject::connect(ui.doubleSpinBoxColorSoftness, SIGNAL(valueChanged(double)), this, SLOT(ChangeColorModelSigmaScale(double)));

	dasp_params_->base_radius = 0.001f * ui.doubleSpinBoxRadius->value();
	dasp_params_->iterations = ui.spinBoxIterations->value();
	dasp_params_->weight_color = ui.doubleSpinBoxWeightColor->value();
	dasp_params_->weight_spatial = ui.doubleSpinBoxWeightSpatial->value();
	dasp_params_->weight_depth = ui.doubleSpinBoxWeightDepth->value();
	dasp_params_->weight_normal = ui.doubleSpinBoxWeightNormal->value();
	dasp_params_->coverage = ui.doubleSpinBoxCoverage->value();

}

WdgtSuperpixelParameters::~WdgtSuperpixelParameters()
{

}

void WdgtSuperpixelParameters::OnSuperSeedType(const QString& txt)
{
	if(txt == "Depth Invariant") {
		dasp_params_->seed_mode = dasp::SeedModes::EquiDistant;
	}
	else if(txt == "Mipmap") {
		dasp_params_->seed_mode = dasp::SeedModes::DepthMipmap;
	}
	else if(txt == "Blue Noise") {
		dasp_params_->seed_mode = dasp::SeedModes::DepthBlueNoise;
	}
	else {
		throw 0;
	}
}

void WdgtSuperpixelParameters::ChangeSuperpixelRadius(double val)
{
	dasp_params_->base_radius = 0.001f * val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelIterations(int val)
{
	dasp_params_->iterations = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightColor(double val)
{
	dasp_params_->weight_color = val;
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

void WdgtSuperpixelParameters::OnColorModelTrain()
{
	if(on_train_) {
		on_train_();
	}
}

void WdgtSuperpixelParameters::ChangeColorModelSigmaScale(double val)
{
	if(on_change_cm_sigma_scale_) {
		on_change_cm_sigma_scale_(val);
	}
}
