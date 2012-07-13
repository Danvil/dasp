#ifndef WDGTBENCHMARK_H
#define WDGTBENCHMARK_H

#include <QtGui/QWidget>
#include <QtGui/QLabel>
#include <QtCore/QTimer>
#include "ui_WdgtBenchmark.h"
#include <Danvil/Tools/Benchmark.h>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <map>
#include <queue>

class WdgtBenchmark : public QWidget
{
    Q_OBJECT

public:
    WdgtBenchmark(QWidget *parent = 0);
    ~WdgtBenchmark();

    void update(const std::string& name, const Danvil::Benchmark::Data& data);

public Q_SLOTS:
	void onTick();

private:
    struct Line
    {
    	boost::shared_ptr<QLabel> l_name;
    	boost::shared_ptr<QLabel> l_mean;
    	boost::shared_ptr<QLabel> l_last;
    };

    Line addLine(bool is_header=false);

private:
    Ui::WdgtBenchmarkClass ui;

    QTimer timer_;

    Line header_;
    std::map<std::string,Line> lines_;

    unsigned int next_line_;

    std::queue<std::pair<std::string,Danvil::Benchmark::Data>> queue_;

    boost::mutex queue_mutex_;
};

#endif // WDGTBENCHMARK_H
