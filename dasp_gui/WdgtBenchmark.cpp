#include "WdgtBenchmark.h"
#include <boost/interprocess/sync/scoped_lock.hpp>

WdgtBenchmark::WdgtBenchmark(QWidget *parent)
    : QWidget(parent)
{
	ui.setupUi(this);
	next_line_ = 0;
	connect(&timer_, SIGNAL(timeout()), this, SLOT(onTick()));
	timer_.setInterval(10);
	timer_.start();

	header_ = addLine(true);
	header_.l_name->setText("Name");
	header_.l_mean->setText("Mean [ms]");
	header_.l_last->setText("Last [ms]");
}

WdgtBenchmark::~WdgtBenchmark()
{

}

void WdgtBenchmark::update(const std::string& name, const Danvil::Benchmark::Data& data)
{
	boost::interprocess::scoped_lock<boost::mutex> lock(queue_mutex_);
	queue_.push({name,data});
}

void WdgtBenchmark::onTick()
{
	boost::interprocess::scoped_lock<boost::mutex> lock(queue_mutex_);
	while(!queue_.empty()) {
		const std::string& name = queue_.front().first;
		const Danvil::Benchmark::Data& data = queue_.front().second;
		auto x = lines_.find(name);
		if(x == lines_.end()) {
			lines_[name] = addLine();
			x = lines_.find(name);
		}
		x->second.l_name->setText(QString::fromStdString(name));
		x->second.l_mean->setText(QString("%1").arg(data.decay_mean,8,'f',3,' '));
		x->second.l_last->setText(QString("%1").arg(data.last,8,'f',3,' '));
		queue_.pop();
	}
}

WdgtBenchmark::Line WdgtBenchmark::addLine(bool is_header)
{
	static QFont fonts[2] = {
			QFont("mono",10,QFont::Normal),
			QFont("mono",10,QFont::Bold)
	};
	Line line;
	line.l_name.reset(new QLabel());
	line.l_name->setFont(fonts[is_header]);
	line.l_mean.reset(new QLabel());
	line.l_mean->setFont(fonts[is_header]);
	line.l_last.reset(new QLabel());
	line.l_last->setFont(fonts[is_header]);
	ui.gridLayout->addWidget(line.l_name.get(),next_line_, 0);
	ui.gridLayout->addWidget(line.l_mean.get(),next_line_, 1);
	ui.gridLayout->addWidget(line.l_last.get(),next_line_, 2);
	next_line_++;
	return line;
}
