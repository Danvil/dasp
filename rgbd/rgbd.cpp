#include "rgbd.hpp"
#include <Slimage/IO.hpp>
#include <Slimage/Gui.hpp>
#include <common/color.hpp>
#include <boost/filesystem.hpp>
#ifdef DASP_HAS_OPENNI
	#include "KinectGrabber.h"
#endif
#include <iostream>
#include <stdexcept>
#include <random>
#include <iterator>
#include <fstream>
#include <snappy.h>

constexpr int WIDTH = 640;
constexpr int HEIGHT = 480;

std::mt19937 random_engine;

slimage::Pixel3ub color(int r, int g, int b)
{
	return {{
		static_cast<unsigned char>(std::min(255, std::max(0, r))),
		static_cast<unsigned char>(std::min(255, std::max(0, g))),
		static_cast<unsigned char>(std::min(255, std::max(0, b)))}};
}

template<typename Dist>
slimage::Pixel3ub color_with_noise(int r, int g, int b, Dist d)
{
	return color(
		r + d(random_engine),
		g + d(random_engine),
		b + d(random_engine));
}

template<bool UseNoise=true>
class RgbdStreamTestUniform : public RgbdStream
{
public:
	RgbdStreamTestUniform() {
		data_.color = slimage::Image3ub(WIDTH, HEIGHT);
		data_.depth = slimage::Image1ui16(WIDTH, HEIGHT);
		constexpr int NC = UseNoise ? 15 : 0;
		constexpr int ND = UseNoise ? 50 : 0;
		std::uniform_int_distribution<int> dc(-NC,+NC);
		for(int i=0; i<HEIGHT; i++) {
			for(int j=0; j<WIDTH; j++) {
				int dd = static_cast<int>(
					ND*std::sin(30.0f*static_cast<float>(i)/static_cast<float>(HEIGHT))
					*std::sin(30.0f*static_cast<float>(j)/static_cast<float>(WIDTH))
				);
				data_.depth(j,i) = 1200 + dd;
				data_.color(j,i) = color_with_noise(32,128,128, dc);
			}
		}
	}
	bool grab() { return true; }
	Rgbd get() { return data_; }
private:
	Rgbd data_;
};

class RgbdStreamTestParaboloid : public RgbdStream
{
public:
	RgbdStreamTestParaboloid() {
		data_.color = slimage::Image3ub(WIDTH, HEIGHT, {{0,128,128}});
		data_.depth = slimage::Image1ui16(WIDTH, HEIGHT);
		for(int i=0; i<HEIGHT; i++) {
			int di = i - HEIGHT/2;
			for(int j=0; j<WIDTH; j++) {
				int dj = j - WIDTH/2;
				data_.depth(j,i) = static_cast<uint16_t>(600 + (di*di + dj*dj)/80);
			}
		}
	}
	bool grab() { return true; }
	Rgbd get() { return data_; }
private:
	Rgbd data_;
};

template<bool UseNoise=true>
class RgbdStreamTestSphere : public RgbdStream
{
public:
	RgbdStreamTestSphere() : time_(0) {}
	bool grab() { return true; }
	Rgbd get() {
		slimage::Image3ub color(WIDTH, HEIGHT);
		slimage::Image1ui16 depth(WIDTH, HEIGHT);
		constexpr int NC = UseNoise ? 15 : 0;
		constexpr int ND = UseNoise ? 50 : 0;
		std::uniform_int_distribution<int> dc(-NC,+NC);
		for(int i=0; i<HEIGHT; i++) {
			for(int j=0; j<WIDTH; j++) {
				float phi = 2.0f*3.1415f*static_cast<float>(time_)/100;
				const float r_move = 70.0f;
				int cx = WIDTH/2 + r_move*std::cos(phi);
				int cy = HEIGHT/2 + r_move*std::sin(phi);
				int dj = j - cx;
				int di = i - cy;
				int r = std::sqrt(di*di + dj*dj);
				int dd = static_cast<int>(
					ND*std::sin(30.0f*static_cast<float>(i)/static_cast<float>(HEIGHT))
					*std::sin(30.0f*static_cast<float>(j)/static_cast<float>(WIDTH))
				);
				if(r < 80) {
					depth(j,i) = 800 + dd;
					color(j,i) = color_with_noise(192,32,32, dc);
				}
				else {
					depth(j,i) = 1200 + dd;
					color(j,i) = color_with_noise(32,128,128, dc);
				}
			}
		}
		time_++;
		return {color, depth};
	}
private:
	int time_;
};

class RgbdStreamStatic : public RgbdStream
{
public:
	RgbdStreamStatic(const std::string& fn) {
		data_.color = slimage::Load3ub(fn + "_color.png");
		data_.depth = slimage::Load1ui16(fn + "_depth.pgm");
		if(data_.color.width() != data_.depth.width() || data_.color.height() != data_.depth.height()) {
			std::cerr << "WARNING: Size of color and depth image do not match!" << std::endl;
		}
	}
	bool grab() { return true; }
	Rgbd get() { return data_; }
private:
	Rgbd data_;
};

class RgbdStreamImages : public RandomAccessRgbdStream
{
public:
	RgbdStreamImages(const std::string& fn)
	:base_dir_(fn) {
		// prepare images
		img_depth_ = slimage::Image1ui16(640, 480, slimage::Pixel1ui16{1000});
		img_color_ = slimage::Image3ub(640, 480, {{0,0,0}});
		// create index
		create_index();
		cur_frame_ = 0;
	}
	~RgbdStreamImages() {
	}
	unsigned int numFrames() {
		return num_frames_;
	}
	unsigned int tell() {
		return cur_frame_;
	}
	void seek(unsigned int frame) {
		cur_frame_ = frame;
	}
	bool grab() {
		if(cur_frame_ < num_frames_) {
			// std::cout << "Current frame " << cur_frame_ << std::endl;
			img_depth_ = slimage::Load1ui16(base_dir_ + "/" + index_[cur_frame_].first);
			img_color_ = slimage::Load3ub(base_dir_ + "/" + index_[cur_frame_].second);
			cur_frame_ ++;
			return true;
		}
		else {
			cur_frame_ = num_frames_;
			return false;
		}
	}
	Rgbd get() {
		return {
			img_color_.clone(),
			img_depth_.clone()
		};
	}
private:
	void create_index() {
		namespace bfs = boost::filesystem;
		bfs::path dir_path(base_dir_);
		bfs::directory_iterator it_end;
		bfs::directory_iterator it(dir_path);
		std::vector<std::string> fns_depth;
		std::vector<std::string> fns_color; 
		for(; it!=it_end; ++it) {
			if(bfs::is_directory(it->status())) {
				continue;
			}
			// pgm -> depth
			if(it->path().extension() == ".pgm") {
				fns_depth.push_back(it->path().filename().string());
			}
			// ppm -> color
			if(it->path().extension() == ".ppm") {
				fns_color.push_back(it->path().filename().string());
			}
		}
		if(fns_depth.size() != fns_color.size()) {
			std::cerr << "Number of color and depht images do not match" << std::endl;
		}
		num_frames_ = fns_depth.size();
		std::sort(fns_depth.begin(), fns_depth.end());
		std::sort(fns_color.begin(), fns_color.end());
		for(unsigned int i=0; i<num_frames_; i++) {
			index_.push_back({fns_depth[i], fns_color[i]});
		}
		std::cout << "Images stream with " << num_frames_ << " frames" << std::endl;
	}
private:
	std::vector<std::pair<std::string,std::string>> index_;
	std::string base_dir_;
	slimage::Image3ub img_color_;
	slimage::Image1ui16 img_depth_;
	unsigned int num_frames_;
	unsigned int cur_frame_;
};

class RgbdStreamFreenectRecord : public RandomAccessRgbdStream
{
public:
	RgbdStreamFreenectRecord(const std::string& fn)
	:base_dir_(fn) {
		// prepare images
		img_depth_ = slimage::Image1ui16(640, 480, slimage::Pixel1ui16{1000});
		img_color_ = slimage::Image3ub(640, 480, {{0,0,0}});
		// read index file
		std::ifstream f_index(base_dir_ + "/INDEX.txt");
		if(!f_index.is_open()) {
			throw std::runtime_error("Could not open index file!");
		}
		std::vector<std::string> raw_index;
		while(f_index) {
			std::string q;
			f_index >> q;
			raw_index.push_back(q);
		}
		parse_index(raw_index);
		cur_frame_ = 0;
	}
	~RgbdStreamFreenectRecord() {
	}
	unsigned int numFrames() {
		return num_frames_;
	}
	unsigned int tell() {
		return cur_frame_;
	}
	void seek(unsigned int frame) {
		cur_frame_ = frame;
	}
	bool grab() {
		cur_frame_ ++;
		if(cur_frame_ < num_frames_) {
			std::cout << "Current frame " << cur_frame_ << std::endl;
			unsnap_depth(base_dir_ + "/" + index_[cur_frame_].first);
			unsnap_color(base_dir_ + "/" + index_[cur_frame_].second);
			return true;
		}
		else {
			cur_frame_ = num_frames_;
			return false;
		}
	}
	Rgbd get() {
		return {
			img_color_.clone(),
			img_depth_.clone()
		};
	}
private:
	void parse_index(const std::vector<std::string>& raw_index) {
		int i=0;
		while(i < raw_index.size()) {
			// wait for depth frame
			int i_depth = -1;
			while(i < raw_index.size()) {
				std::string q = raw_index[i];
				if(q[0] == 'd') {
					i_depth = i;
					break;
				}
				i++;
			}
			// wait for color frame
			int i_color = -1;
			while(i < raw_index.size()) {
				std::string q = raw_index[i];
				if(q[0] == 'r') {
					i_color = i;
					break;
				}
				i++;
			}
			// add to index
			if(i_depth == -1 || i_color == -1) {
				break;
			}
			index_.push_back({raw_index[i_depth], raw_index[i_color]});
		}
		num_frames_ = index_.size();
		std::cout << "freenect stream with " << num_frames_ << " frames" << std::endl;
	}
	std::vector<char> unsnap(const std::string& fn) {
		// read data
		std::ifstream ifs(fn);
		if(!ifs.is_open()) {
			throw std::runtime_error("Could not open file!");
		}
		ifs.seekg(0, std::ios::end);
		std::size_t compressed_length = ifs.tellg();
		std::cout << "compressed_length "  << compressed_length << std::endl;
		ifs.seekg(0, std::ios::beg);
		std::vector<char> compressed(compressed_length);
		ifs.read(compressed.data(), compressed_length);
		// unsnap
		std::size_t uncompressed_length = 0;
		snappy::GetUncompressedLength(compressed.data(), compressed.size(), &uncompressed_length);
		std::cout << "uncompressed_length "  << uncompressed_length << std::endl;
		std::vector<char> uncompressed(uncompressed_length);
		snappy::RawUncompress(compressed.data(), compressed.size(), uncompressed.data());
		return uncompressed;
	}
	void unsnap_depth(const std::string& fn) {
		std::cout << "unsnap depth" << std::endl;
		std::vector<char> v = unsnap(fn);
		std::copy(v.begin(), v.end(), (char*)(img_depth_.begin().pointer()));
		// slimage::gui::Show("depth", common::GreyDepth(img_depth_,500,3000), 10);
		std::cout << "depth ready" << std::endl;
	}
	void unsnap_color(const std::string& fn) {
		std::cout << "unsnap color" << std::endl;
		std::vector<char> v = unsnap(fn);
		std::copy(v.begin(), v.end(), img_color_.begin().pointer());
		// slimage::gui::Show("color", img_color_, 10);
		std::cout << "color ready" << std::endl;
	}
private:
	std::vector<std::pair<std::string,std::string>> index_;
	std::string base_dir_;
	slimage::Image3ub img_color_;
	slimage::Image1ui16 img_depth_;
	unsigned int num_frames_;
	unsigned int cur_frame_;
};

#ifdef DASP_HAS_OPENNI

class RgbdStreamOni : public RandomAccessRgbdStream
{
public:
	RgbdStreamOni(const std::string& fn) {
		grabber_ = std::make_shared<dasp::KinectGrabber>();
		grabber_->OpenFile(fn);
	}
	~RgbdStreamOni() {
		grabber_->Stop();
	}
	unsigned int numFrames() {
		return grabber_->NumFrames();
	}
	unsigned int tell() {
		return grabber_->TellFrame();
	}
	void seek(unsigned int frame) {
		grabber_->SeekToFrame(frame);
	}
	bool grab() {
		return grabber_->Grab();
	}
	Rgbd get() {
		return {
			grabber_->GetLastColor().clone(),
			grabber_->GetLastDepth().clone()
		};
	}
private:
	std::shared_ptr<dasp::KinectGrabber> grabber_;
};

class RgbdStreamKinectLive : public RgbdStream
{
public:
	RgbdStreamKinectLive(const std::string& fn_config) {
		grabber_ = std::make_shared<dasp::KinectGrabber>();
		grabber_->OpenConfig(fn_config);
	}
	~RgbdStreamKinectLive() {
		grabber_->Stop();
	}
	bool grab() {
		return grabber_->Grab();
	}
	Rgbd get() {
		return {
			grabber_->GetLastColor().clone(),
			grabber_->GetLastDepth().clone()
		};
	}
private:
	std::shared_ptr<dasp::KinectGrabber> grabber_;
};

#endif

std::shared_ptr<RgbdStream> FactorTest(const std::string& arg)
{
	if(arg == "uniform") {
		return std::make_shared<RgbdStreamTestUniform<false>>();
	}
	if(arg == "paraboloid") {
		return std::make_shared<RgbdStreamTestParaboloid>();
	}
	if(arg == "sphere") {
		return std::make_shared<RgbdStreamTestSphere<false>>();
	}
	std::cerr << "ERROR: Invalid rgbd test stream arg='" << arg << "'!" << std::endl;
	throw 0;
}

std::shared_ptr<RgbdStream> FactorStatic(const std::string& fn)
{
	return std::make_shared<RgbdStreamStatic>(fn);
}

std::shared_ptr<RandomAccessRgbdStream> FactorImages(const std::string& fn)
{
	return std::make_shared<RgbdStreamImages>(fn);
}

std::shared_ptr<RandomAccessRgbdStream> FactorFreenectRecord(const std::string& fn)
{
	return std::make_shared<RgbdStreamFreenectRecord>(fn);
}

std::shared_ptr<RandomAccessRgbdStream> FactorOni(const std::string& fn)
{
#ifdef DASP_HAS_OPENNI
	return std::make_shared<RgbdStreamOni>(fn);
#else
	std::cerr << "ERROR: Library configured without OpenNI! Enabel DASP_HAS_OPENNI in CMake." << std::endl;
	throw 0;
#endif
}

std::shared_ptr<RgbdStream> FactorKinectLive(const std::string& fn_config)
{
#ifdef DASP_HAS_OPENNI
	return std::make_shared<RgbdStreamKinectLive>(fn_config);
#else
	std::cerr << "ERROR: Library configured without OpenNI! Enabel DASP_HAS_OPENNI in CMake." << std::endl;
	throw 0;
#endif
}

std::shared_ptr<RgbdStream> FactorStream(const std::string& mode, const std::string& arg)
{
	if(mode == "test") {
		return FactorTest(arg);
	}
	if(mode == "static") {
		return FactorStatic(arg);
	}
	if(mode == "images") {
		return FactorImages(arg);
	}
	if(mode == "freenect") {
		return FactorFreenectRecord(arg);
	}
	if(mode == "oni") {
		return FactorOni(arg);
	}
	if(mode == "live") {
		return FactorKinectLive(arg);
	}
	// default
	std::cerr << "ERROR: Invalid mode='" << mode << "' and arg='" << arg << "'!" << std::endl;
	throw 0;
}
