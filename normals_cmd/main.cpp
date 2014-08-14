#include <dasp/Tools.hpp>
#include <slimage/opencv.hpp>
#include <slimage/io.hpp>
#include <slimage/image.hpp>
#include <Danvil/Tools/MoreMath.h>
#include <boost/program_options.hpp>
#include <iostream>

inline Eigen::Vector3f GradientToNormal(const Eigen::Vector3f& position, const Eigen::Vector2f& g)
{
	const float gx = g.x();
	const float gy = g.y();
	const float scl = Danvil::MoreMath::FastInverseSqrt(gx*gx + gy*gy + 1.0f);
	Eigen::Vector3f normal(scl*gx, scl*gy, -scl);
	// force normal to look towards the camera
	// check if point to camera direction and normal are within 90 deg
	// enforce: normal * (cam_pos - pos) > 0
	// do not need to normalize (cam_pos - pos) as only sign is considered
	const float q = normal.dot(-position);
	if(q < 0) {
		return normal * -1.0f;
	}
	else if(q == 0) {
		// this should not happen ...
		return Eigen::Vector3f(0,0,-1);
	}
	else {
		return normal;
	}
}

std::vector<Eigen::Vector3f> ComputeNormals(const slimage::Image1ui16& depth, const dasp::Camera& cam)
{
	constexpr float NORMAL_RADIUS = 0.025f;
	int width = depth.width();
	int height = depth.height();
	std::vector<Eigen::Vector3f> normals(width*height, Eigen::Vector3f{0,0,-1});
	for(int y=0; y<height; y++) {
		for(int x=0; x<width; x++) {
			const int i = x + y*width;
			const uint16_t depth_i16 = depth[i];
			if(depth_i16 == 0) {
				continue;
			}
			const float z_over_f = cam.convertKinectToMeter(depth_i16) / cam.focal;
			const float cluster_radius_px = NORMAL_RADIUS / z_over_f;
			Eigen::Vector2f gradient = LocalDepthGradient(depth, x, y, z_over_f, 0.5f*cluster_radius_px, cam);
			Eigen::Vector3f position = cam.unprojectImpl(static_cast<float>(x), static_cast<float>(y), z_over_f);
			Eigen::Vector3f normal = GradientToNormal(position, gradient);
			normals[i] = normal;
		}
	}
	return normals;
}

inline slimage::Pixel3ub ColorizeNormal(const Eigen::Vector3f& n)
{
	return slimage::Pixel3ub{
		static_cast<unsigned char>(255.0f * 0.5f * (1.0f + n.x())),
		static_cast<unsigned char>(255.0f * 0.5f * (1.0f + n.y())),
		static_cast<unsigned char>(255.0f * 0.5f * (1.0f + n.z()))
	};
}

slimage::Image3ub ColorizeNormals(const std::vector<Eigen::Vector3f>& normals, int width, int height)
{
	slimage::Image3ub img(width, height);
	for(int i=0; i<img.size(); i++) {
		img[i] = ColorizeNormal(normals[i]);
	}
	return img;
}

int main(int argc, char** argv)
{
	std::string p_img = "";
	std::string p_out = "out";
	bool p_verbose = false;
	dasp::Camera p_cam{320.0f, 240.0f, 540.0f, 0.001f};

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("img", po::value(&p_img), "path to depth image (must be a pgm file)")
		("out", po::value(&p_out), "path to result image with color encoded normals")
		("verbose", po::value(&p_verbose), "verbose")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	// loading depth
	if(p_verbose) std::cout << "Reading DEPTH: '" << p_img << "'" << std::endl;
	slimage::Image1ui16 img_depth =  slimage::Load1ui16(p_img);

	// computing normals
	std::vector<Eigen::Vector3f> normals = ComputeNormals(img_depth, p_cam);

	// colorize
	slimage::Image3ub img_normals = ColorizeNormals(normals, img_depth.width(), img_depth.height());
	slimage::Save(p_out, img_normals);

	return 0;
}
