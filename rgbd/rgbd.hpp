#ifndef RGBD_HPP
#define RGBD_HPP

#include <slimage/image.hpp>
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

	Rgbd grabAndGet(bool* x=0) {
		bool v = grab();
		if(x) *x = v;
		return get();
	}
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

std::shared_ptr<RgbdStream> FactorTest(const std::string& arg);

std::shared_ptr<RgbdStream> FactorStatic(const std::string& fn);

std::shared_ptr<RandomAccessRgbdStream> FactorImages(const std::string& fn);

std::shared_ptr<RandomAccessRgbdStream> FactorFreenectRecord(const std::string& fn);

std::shared_ptr<RandomAccessRgbdStream> FactorOni(const std::string& fn);

std::shared_ptr<RgbdStream> FactorKinectLive(const std::string& fn_config);

std::shared_ptr<RgbdStream> FactorStream(const std::string& mode, const std::string& arg);

#endif
