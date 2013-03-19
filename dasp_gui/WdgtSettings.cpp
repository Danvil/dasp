#include "WdgtSettings.h"

WdgtSettings::WdgtSettings(const boost::shared_ptr<DaspProcessing>& dasp_processing, QWidget *parent)
    : QWidget(parent)
{
	dasp_processing_ = dasp_processing;

	ui.setupUi(this);

	ui.comboBoxPlotGraphWeights->addItem("White", 0);
	ui.comboBoxPlotGraphWeights->addItem("Black", 1);
	ui.comboBoxPlotGraphWeights->addItem("DASP metric", 2);
	ui.comboBoxPlotGraphWeights->addItem("Spectral metric", 3);
	ui.comboBoxPlotGraphWeights->addItem("Spectral result", 4);
	ui.comboBoxPlotGraphWeights->setCurrentIndex(2);

	ui.comboBoxPlotPointsColor->addItem("Color", dasp::plots::Color);
	ui.comboBoxPlotPointsColor->addItem("Depth", dasp::plots::Depth);
	ui.comboBoxPlotPointsColor->addItem("Valid", dasp::plots::Valid);
	ui.comboBoxPlotPointsColor->addItem("Gradient", dasp::plots::Gradient);
	ui.comboBoxPlotPointsColor->addItem("ClusterRadius", dasp::plots::ClusterRadius);
	ui.comboBoxPlotPointsColor->addItem("Density", dasp::plots::Density);
	ui.comboBoxPlotPointsColor->setCurrentIndex(0);

	ui.comboBoxPlotClusterColor->addItem("UniBlack", dasp::plots::UniBlack);
	ui.comboBoxPlotClusterColor->addItem("UniWhite", dasp::plots::UniWhite);
	ui.comboBoxPlotClusterColor->addItem("Color", dasp::plots::Color);
	ui.comboBoxPlotClusterColor->addItem("Depth", dasp::plots::Depth);
	ui.comboBoxPlotClusterColor->addItem("Gradient", dasp::plots::Gradient);
	ui.comboBoxPlotClusterColor->addItem("Thickness", dasp::plots::Thickness);
	ui.comboBoxPlotClusterColor->addItem("Eccentricity", dasp::plots::Eccentricity);
	ui.comboBoxPlotClusterColor->addItem("AreaQuotient", dasp::plots::AreaQuotient);
	ui.comboBoxPlotClusterColor->addItem("CoverageError", dasp::plots::CoverageError);
	ui.comboBoxPlotClusterColor->setCurrentIndex(2);

	ui.comboBoxPlotClusterMode->addItem("ClusterCenter", dasp::plots::ClusterCenter);
	ui.comboBoxPlotClusterMode->addItem("ClusterPoints", dasp::plots::ClusterPoints);
	ui.comboBoxPlotClusterMode->addItem("ClusterEllipses", dasp::plots::ClusterEllipses);
	ui.comboBoxPlotClusterMode->addItem("ClusterEllipsesFilled", dasp::plots::ClusterEllipsesFilled);
	ui.comboBoxPlotClusterMode->setCurrentIndex(1);

	QObject::connect(ui.doubleSpinBoxDaspSegmentThreshold, SIGNAL(valueChanged(double)), this, SLOT(ChangeDaspSegmentThreshold(double)));
	QObject::connect(ui.checkBoxPlotPoints, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotPoints(int)));
	QObject::connect(ui.comboBoxPlotPointsColor, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangePlotPointsColor(int)));
	QObject::connect(ui.checkBoxPlotClusters, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotClusters(int)));
	QObject::connect(ui.comboBoxPlotClusterMode, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangePlotClusterMode(int)));
	QObject::connect(ui.comboBoxPlotClusterColor, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangePlotClusterColor(int)));
	QObject::connect(ui.checkBoxPlotBorders, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotBorders(int)));
	QObject::connect(ui.checkBoxPlotGraphSpatialCut, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotGraphSpatialCut(int)));
	QObject::connect(ui.checkBoxPlotGraph, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotGraph(int)));
	QObject::connect(ui.comboBoxPlotGraphWeights, SIGNAL(currentIndexChanged(int)), this, SLOT(ChangePlotGraphWeights(int)));
	QObject::connect(ui.checkBoxPlotDensity, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotDensity(int)));
	QObject::connect(ui.checkBoxPlotSegments, SIGNAL(stateChanged(int)), this, SLOT(ChangePlotSegments(int)));

	dasp_processing_->dasp_params->segment_threshold = ui.doubleSpinBoxDaspSegmentThreshold->value();
	dasp_processing_->show_points_ = ui.checkBoxPlotPoints->isChecked();
	dasp_processing_->point_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotPointsColor->itemData(ui.comboBoxPlotPointsColor->currentIndex()).toInt());
	dasp_processing_->show_clusters_ = ui.checkBoxPlotClusters->isChecked();;
	dasp_processing_->cluster_mode_ = (dasp::plots::ClusterMode)(ui.comboBoxPlotClusterMode->itemData(ui.comboBoxPlotClusterMode->currentIndex()).toInt());
	dasp_processing_->cluster_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotClusterColor->itemData(ui.comboBoxPlotClusterColor->currentIndex()).toInt());
	dasp_processing_->show_cluster_borders_ = ui.checkBoxPlotBorders->isChecked();
	dasp_processing_->graph_cut_spatial_ = ui.checkBoxPlotGraphSpatialCut->isChecked();
	dasp_processing_->show_graph_ = ui.checkBoxPlotGraph->isChecked();
	dasp_processing_->show_graph_weights_ = ui.comboBoxPlotGraphWeights->itemData(ui.comboBoxPlotGraphWeights->currentIndex()).toInt();
	dasp_processing_->plot_density_ = ui.checkBoxPlotDensity->isChecked();
	dasp_processing_->plot_segments_ = ui.checkBoxPlotSegments->isChecked();

}

WdgtSettings::~WdgtSettings()
{

}

void WdgtSettings::ChangeDaspSegmentThreshold(double val)
{
	dasp_processing_->dasp_params->segment_threshold = val;
	*reload = true;
}

void WdgtSettings::ChangePlotPoints(int state)
{
	dasp_processing_->show_points_ = state;
	*reload = true;
}

void WdgtSettings::ChangePlotPointsColor(int selection)
{
	dasp_processing_->point_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotPointsColor->itemData(selection).toInt());
	*reload = true;
}

void WdgtSettings::ChangePlotClusters(int state)
{
	dasp_processing_->show_clusters_ = state;
	*reload = true;
}

void WdgtSettings::ChangePlotClusterMode(int selection)
{
	dasp_processing_->cluster_mode_ = (dasp::plots::ClusterMode)(ui.comboBoxPlotClusterMode->itemData(selection).toInt());
	*reload = true;
}

void WdgtSettings::ChangePlotClusterColor(int selection)
{
	dasp_processing_->cluster_color_mode_ = (dasp::plots::ColorMode)(ui.comboBoxPlotClusterColor->itemData(selection).toInt());
	*reload = true;
}

void WdgtSettings::ChangePlotBorders(int state)
{
	dasp_processing_->show_cluster_borders_ = state;
	*reload = true;
}

void WdgtSettings::ChangePlotGraphSpatialCut(int state)
{
	dasp_processing_->graph_cut_spatial_ = state;
	*reload = true;
}

void WdgtSettings::ChangePlotGraph(int state)
{
	dasp_processing_->show_graph_ = state;
	*reload = true;
}

void WdgtSettings::ChangePlotGraphWeights(int selection)
{
	dasp_processing_->show_graph_weights_ = ui.comboBoxPlotGraphWeights->itemData(selection).toInt();
	*reload = true;
}

void WdgtSettings::ChangePlotDensity(int state)
{
	dasp_processing_->plot_density_ = state;
	*reload = true;
}

void WdgtSettings::ChangePlotSegments(int state)
{
	dasp_processing_->plot_segments_ = state;
	*reload = true;
}
