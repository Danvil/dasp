#include "WdgtSuperpixelParameters.h"

WdgtSuperpixelParameters::WdgtSuperpixelParameters(const boost::shared_ptr<dasp::DaspTracker>& dasp_tracker, QWidget *parent)
    : QWidget(parent)
{
	dasp_tracker_ = dasp_tracker;

	ui.setupUi(this);

	ui.comboBoxSeedType->addItem("EquiDistant", dasp::SeedModes::EquiDistant);
	ui.comboBoxSeedType->addItem("DepthMipmap", dasp::SeedModes::DepthMipmap);
	ui.comboBoxSeedType->addItem("DepthBlueNoise", dasp::SeedModes::DepthBlueNoise);
	ui.comboBoxSeedType->addItem("Delta", dasp::SeedModes::Delta);
	ui.comboBoxSeedType->setCurrentIndex(1);

	ui.comboBoxPlotPointsColor->addItem("Color", dasp::plots::Color);
	ui.comboBoxPlotPointsColor->addItem("Depth", dasp::plots::Depth);
	ui.comboBoxPlotPointsColor->addItem("Gradient", dasp::plots::Gradient);
	ui.comboBoxPlotPointsColor->addItem("Color", dasp::plots::Color);
	ui.comboBoxPlotPointsColor->setCurrentIndex(0);

	ui.comboBoxPlotClusterColor->addItem("UniBlack", dasp::plots::UniBlack);
	ui.comboBoxPlotClusterColor->addItem("UniWhite", dasp::plots::UniWhite);
	ui.comboBoxPlotClusterColor->addItem("Color", dasp::plots::Color);
	ui.comboBoxPlotClusterColor->addItem("Depth", dasp::plots::Depth);
	ui.comboBoxPlotClusterColor->addItem("Gradient", dasp::plots::Gradient);
	ui.comboBoxPlotClusterColor->addItem("Color", dasp::plots::Color);
	ui.comboBoxPlotClusterColor->setCurrentIndex(2);

	ui.comboBoxPlotClusterMode->addItem("ClusterCenter", dasp::plots::ClusterCenter);
	ui.comboBoxPlotClusterMode->addItem("ClusterPoints", dasp::plots::ClusterPoints);
	ui.comboBoxPlotClusterMode->addItem("ClusterEllipses", dasp::plots::ClusterEllipses);
	ui.comboBoxPlotClusterMode->addItem("ClusterEllipsesFilled", dasp::plots::ClusterEllipsesFilled);
	ui.comboBoxPlotClusterMode->setCurrentIndex(1);

	QObject::connect(ui.comboBoxSeedType, SIGNAL(currentIndexChanged(int)), this, SLOT(OnSuperSeedType(int)));
	QObject::connect(ui.checkBoxGradientAdaptive, SIGNAL(stateChanged(int)), this, SLOT(ChangeSuperUseGradientDensity(int)));
	QObject::connect(ui.doubleSpinBoxRadius, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelRadius(double)));
	QObject::connect(ui.spinBoxIterations, SIGNAL(valueChanged(int)), this, SLOT(ChangeSuperpixelIterations(int)));
	QObject::connect(ui.doubleSpinBoxWeightColor, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightColor(double)));
	QObject::connect(ui.doubleSpinBoxWeightSpatial, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightSpatial(double)));
	QObject::connect(ui.doubleSpinBoxWeightDepth, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightDepth(double)));
	QObject::connect(ui.doubleSpinBoxWeightNormal, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelWeightNormal(double)));
	QObject::connect(ui.doubleSpinBoxCoverage, SIGNAL(valueChanged(double)), this, SLOT(ChangeSuperpixelCoverage(double)));
	QObject::connect(ui.pushButtonColorModelTrain, SIGNAL(clicked()), this, SLOT(OnColorModelTrain()));
	QObject::connect(ui.doubleSpinBoxColorSoftness, SIGNAL(valueChanged(double)), this, SLOT(ChangeColorModelSigmaScale(double)));

	QObject::connect(ui.checkBoxPlotPoints, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotPoints(int)));
	QObject::connect(ui.comboBoxPlotPointsColor, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangePlotPointsColor(int)));
	QObject::connect(ui.checkBoxPlotClusters, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotClusters(int)));
	QObject::connect(ui.comboBoxPlotClusterMode, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangePlotClusterMode(int)));
	QObject::connect(ui.comboBoxPlotClusterColor, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangePlotClusterColor(int)));
	QObject::connect(ui.checkBoxPlotBorders, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotBorders(int)));
	QObject::connect(ui.checkBoxPlotGraph, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotGraph(int)));

	dasp_tracker_->dasp_params->gradient_adaptive_density = ui.checkBoxGradientAdaptive->isChecked();
	dasp_tracker_->dasp_params->base_radius = 0.001f * ui.doubleSpinBoxRadius->value();
	dasp_tracker_->dasp_params->iterations = ui.spinBoxIterations->value();
	dasp_tracker_->dasp_params->weight_color = ui.doubleSpinBoxWeightColor->value();
	dasp_tracker_->dasp_params->weight_spatial = ui.doubleSpinBoxWeightSpatial->value();
	dasp_tracker_->dasp_params->weight_depth = ui.doubleSpinBoxWeightDepth->value();
	dasp_tracker_->dasp_params->weight_normal = ui.doubleSpinBoxWeightNormal->value();
	dasp_tracker_->dasp_params->coverage = ui.doubleSpinBoxCoverage->value();

	dasp_tracker_->show_points_ = ui.checkBoxPlotPoints->isChecked();;
	dasp_tracker_->point_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotPointsColor->itemData(ui.comboBoxPlotPointsColor->currentIndex()).toInt());
	dasp_tracker_->show_clusters_ = ui.checkBoxPlotClusters->isChecked();;
	dasp_tracker_->cluster_mode_ = (dasp::plots::ClusterMode)(ui.comboBoxPlotClusterMode->itemData(ui.comboBoxPlotClusterMode->currentIndex()).toInt());
	dasp_tracker_->cluster_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotClusterColor->itemData(ui.comboBoxPlotClusterColor->currentIndex()).toInt());
	dasp_tracker_->show_cluster_borders_ = ui.checkBoxPlotBorders->isChecked();
	dasp_tracker_->show_graph_ = ui.checkBoxPlotGraph->isChecked();

}

WdgtSuperpixelParameters::~WdgtSuperpixelParameters()
{

}

void WdgtSuperpixelParameters::OnSuperSeedType(int selection)
{
	dasp_tracker_->dasp_params->seed_mode = (dasp::SeedMode)(ui.comboBoxSeedType->itemData(selection).toInt());
}

void WdgtSuperpixelParameters::ChangeSuperUseGradientDensity(int state)
{
	dasp_tracker_->dasp_params->gradient_adaptive_density = state;
}

void WdgtSuperpixelParameters::ChangeSuperpixelRadius(double val)
{
	dasp_tracker_->dasp_params->base_radius = 0.001f * val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelIterations(int val)
{
	dasp_tracker_->dasp_params->iterations = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightColor(double val)
{
	dasp_tracker_->dasp_params->weight_color = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightSpatial(double val)
{
	dasp_tracker_->dasp_params->weight_spatial = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightDepth(double val)
{
	dasp_tracker_->dasp_params->weight_depth = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelWeightNormal(double val)
{
	dasp_tracker_->dasp_params->weight_normal = val;
}

void WdgtSuperpixelParameters::ChangeSuperpixelCoverage(double val)
{
	dasp_tracker_->dasp_params->coverage = val;
}

void WdgtSuperpixelParameters::OnColorModelTrain()
{
	dasp_tracker_->training_ = true;
}

void WdgtSuperpixelParameters::ChangeColorModelSigmaScale(double val)
{
	dasp_tracker_->color_model_sigma_scale_ = val;
}

void WdgtSuperpixelParameters::ChangePlotPoints(int state)
{
	dasp_tracker_->show_points_ = state;
}

void WdgtSuperpixelParameters::ChangePlotPointsColor(int selection)
{
	dasp_tracker_->point_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotPointsColor->itemData(selection).toInt());
}

void WdgtSuperpixelParameters::ChangePlotClusters(int state)
{
	dasp_tracker_->show_clusters_ = state;
}

void WdgtSuperpixelParameters::ChangePlotClusterMode(int selection)
{
	dasp_tracker_->cluster_mode_ = (dasp::plots::ClusterMode)(ui.comboBoxPlotClusterMode->itemData(selection).toInt());
}

void WdgtSuperpixelParameters::ChangePlotClusterColor(int selection)
{
	dasp_tracker_->cluster_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotClusterColor->itemData(selection).toInt());
}

void WdgtSuperpixelParameters::ChangePlotBorders(int state)
{
	dasp_tracker_->show_cluster_borders_ = state;
}

void WdgtSuperpixelParameters::ChangePlotGraph(int state)
{
	dasp_tracker_->show_graph_ = state;
}
