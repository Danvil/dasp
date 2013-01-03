#ifndef RGBD_HPP
#define RGBD_HPP

#include <Slimage/Slimage.hpp>
#include <memory>

struct Rgbd
{
	slimage::Image3ub color;
	slimage::Image1ui16 depth;
};

class Scenario
{
public:
	virtual ~Scenario() {}
	virtual Rgbd pop() = 0;
};

std::shared_ptr<Scenario> FactorTest(const std::string& tag);

std::shared_ptr<Scenario> FactorStatic(const std::string& fn);

std::shared_ptr<Scenario> FactorOni(const std::string& fn, unsigned int offset=0);

std::shared_ptr<Scenario> FactorLive();

#endif
