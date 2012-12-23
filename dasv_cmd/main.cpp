#include <dasv.hpp>
#include <Slimage/IO.hpp>
#include <kinect/KinectGrabber.h>
#include <string>
#include <iostream>
#include <memory>
#include <Slimage/Gui.hpp>

int main(int argc, char** argv)
{
	using namespace dasv;

//	const std::string ds_path = "/home/david/Documents/DataSets";
	const std::string ds_path = "/home/david/Dokumente/DataSets";

	const int WIDTH = 640;
	const int HEIGHT = 480;

	// std::cout << "Uniform test" << std::endl;
	// slimage::Image3ub color(WIDTH, HEIGHT, {{0,128,128}});
	// slimage::Image1ui16 depth(WIDTH, HEIGHT);
	// for(int i=0; i<HEIGHT; i++) {
	// 	for(int j=0; j<WIDTH; j++) {
	// 		//depth(j,i) = 1000;
	// 		int di = i - HEIGHT/2;
	// 		int dj = j - WIDTH/2;
	// 		depth(j,i) = static_cast<uint16_t>(600 + (di*di + dj*dj)/80);
	// 	}
	// }
	// ContinuousSupervoxels sv;
	// sv.start(WIDTH,HEIGHT);
	// for(int t=0; t<1000; t++) {
	// 	slimage::gui::Show("color", color, 0);
	// 	slimage::gui::Show("depth", depth, 500, 3000, 0);
	// 	sv.step(color, depth);
	// 	slimage::gui::WaitForKeypress();
	// }


	// std::cout << "moving sphere test" << std::endl;
	// ContinuousSupervoxels sv;
	// sv.start(WIDTH,HEIGHT);
	// for(int t=0; t<1000; t++) {
	// 	slimage::Image3ub color(WIDTH, HEIGHT, {{0,128,128}});
	// 	slimage::Image1ui16 depth(WIDTH, HEIGHT);
	// 	for(int i=0; i<HEIGHT; i++) {
	// 		for(int j=0; j<WIDTH; j++) {
	// 			depth(j,i) = 1200;
	// 			float phi = 2.0f*3.1415f*static_cast<float>(t)/100;
	// 			const float r_move = 70.0f;
	// 			int cx = WIDTH/2 + r_move*std::cos(phi);
	// 			int cy = HEIGHT/2 + r_move*std::sin(phi);
	// 			int di = i - cx;
	// 			int dj = j - cy;
	// 			int d = std::sqrt(di*di + dj*dj);
	// 			if(d < 80) {
	// 				depth(j,i) = 800;
	// 				color(j,i) = {{255,0,0}};
	// 			}
	// 		}
	// 	}
	// 	slimage::gui::Show("color", color, 0);
	// 	slimage::gui::Show("depth", depth, 500, 3000, 0);
	// 	sv.step(color, depth);
	// 	slimage::gui::WaitForKeypress();
	// }


	std::cout << "Loading frame data" << std::endl;
	std::string fn_color = ds_path + "/dasp_rgbd_dataset/images/001_color.png";
	std::string fn_depth = ds_path + "/dasp_rgbd_dataset/images/001_depth.pgm";
	slimage::Image3ub color = slimage::Load3ub(fn_color);
//	slimage::gui::Show("color", color, 0);
	slimage::Image1ui16 depth = slimage::Load1ui16(fn_depth);
//	slimage::gui::Show("depth", depth, 500, 3000, 0);	
	ContinuousSupervoxels sv;
	sv.start(WIDTH,HEIGHT);
	for(int t=0; t<200; t++) {
		sv.step(color, depth);
		slimage::gui::WaitForKeypress();
	}


//	std::cout << "Loading ONI data" << std::endl;
// 	std::shared_ptr<dasp::KinectGrabber> grabber = std::make_shared<dasp::KinectGrabber>();
// 	grabber->OpenFile(ds_path + "/2012-10-12 cogwatch dasp/SlowVelocity/C15_c01_slow.oni");
// 	grabber->SeekToFrame(100);
// 	ContinuousSupervoxels sv;
// 	sv.start(WIDTH,HEIGHT);
// 	for(int i=0; i<1000; i++) {
// 		bool ok = grabber->Grab();
// 		if(!ok) {
// 			break;
// 		}
// 		auto color = grabber->GetLastColor().clone();
// 		auto depth = grabber->GetLastDepth().clone();
// 		sv.step(color, depth);
// 	}

	std::vector<Cluster> clusters = sv.getAllClusters();
	std::cout << "Supervoxel count = " << clusters.size() << std::endl;
	DebugWriteClusters("clusters.tsv", clusters);

	std::cout << "Finished." << std::endl;
//	slimage::gui::WaitForKeypress();

	return 1;
}
