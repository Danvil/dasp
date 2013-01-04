#ifndef RGBD_HPP
#define RGBD_HPP

#include <Slimage/Slimage.hpp>
#include <memory>
#include <string>

struct Rgbd
{
	slimage::Image3ub color;
	slimage::Image1ui16 depth;
};

class RgbdStream
{
public:
	virtual ~RgbdStream() {}
	virtual bool grab() = 0;
	virtual Rgbd get() = 0;
};

class RandomAccessRgbdStream
: public RgbdStream
{
public:
	virtual ~RandomAccessRgbdStream() {}
	virtual unsigned int numFrames() = 0;
	virtual unsigned int tell() = 0;
	virtual void seek(unsigned int frame) = 0;
};

std::shared_ptr<RgbdStream> FactorTest(const std::string& tag);

std::shared_ptr<RgbdStream> FactorStatic(const std::string& fn);

std::shared_ptr<RandomAccessRgbdStream> FactorOni(const std::string& fn);

std::shared_ptr<RgbdStream> FactorKinectLive(const std::string& fn_config);

#endif
