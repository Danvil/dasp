#ifndef WDGTABOUT_H
#define WDGTABOUT_H

#include <QtGui/QWidget>
#include <QtGui/QLabel>
#include "ui_WdgtAbout.h"

class WdgtAbout
: public QWidget
{
	Q_OBJECT

public:
	WdgtAbout(QWidget *parent = 0);
	~WdgtAbout();

private:
	Ui::WdgtAboutClass ui;
};

#endif
