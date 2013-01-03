#ifndef DASVVIS_WDGTDASVVIS_H
#define DASVVIS_WDGTDASVVIS_H

#include "ui_WdgtDasvVis.h"
#include <Candy/System/GLSystemQtWindow.h>
#include <Candy.h>
#include <QtGui/QMainWindow>

class WdgtDasvVis : public QMainWindow
{
    Q_OBJECT

public:
	WdgtDasvVis(QWidget *parent = 0);
	~WdgtDasvVis();

private:
 	boost::shared_ptr<Candy::Engine> engine_main_;
	boost::shared_ptr<Candy::Engine> engine_tab_;
	Candy::GLSystemQtWindow* widget_candy_tab_;

private:
    Ui::WdgtDasvVisClass ui;
};

#endif
