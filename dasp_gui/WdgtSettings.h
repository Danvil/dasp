#ifndef DASPGUI_WDGTSETTINGS_H
#define DASPGUI_WDGTSETTINGS_H

#include <QtGui/QWidget>
#include "ui_WdgtSettings.h"
#include "DaspProcessing.h"
#include <functional>

class WdgtSettings : public QWidget
{
    Q_OBJECT

public:
	WdgtSettings(const boost::shared_ptr<DaspProcessing>& dasp_processing, QWidget *parent = 0);
	~WdgtSettings();

	bool* reload;

public Q_SLOTS:
	void ChangeDaspSegmentThreshold(double val);
	void ChangePlotPoints(int state);
	void ChangePlotPointsColor(int selection);
	void ChangePlotClusters(int state);
	void ChangePlotClusterMode(int selection);
	void ChangePlotClusterColor(int selection);
	void ChangePlotBorders(int state);
	void ChangePlotGraphSpatialCut(int state);
	void ChangePlotGraph(int state);
	void ChangePlotGraphWeights(int val);
	void ChangePlotDensity(int state);
	void ChangePlotSegments(int state);

private:
	boost::shared_ptr<DaspProcessing> dasp_processing_;

private:
    Ui::WdgtSettingsClass ui;
};

#endif
