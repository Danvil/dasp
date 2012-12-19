#include <dasv.hpp>
#include <Slimage/IO.hpp>
#include <kinect/KinectGrabber.h>
#include <string>
#include <iostream>
#include <memory>
//#include <Slimage/Gui.hpp>

int main(int argc, char** argv)
{
	using namespace dasv;

// 	std::cout << "Loading data" << std::endl;
// 	std::string fn_color = "/home/david/Documents/DataSets/dasp_rgbd_dataset/images/001_color.png";
// 	std::string fn_depth = "/home/david/Documents/DataSets/dasp_rgbd_dataset/images/001_depth.pgm";
// 	slimage::Image3ub color = slimage::Load3ub(fn_color);
// //	slimage::gui::Show("color", color, 0);
// 	slimage::Image1ui16 depth = slimage::Load1ui16(fn_depth);
// //	slimage::gui::Show("depth", depth, 500, 3000, 0);	
// 	ContinuousSupervoxels sv;
// 	sv.start(640,480);
// 	for(int t=0; t<25; t++) {
// 		sv.step(color, depth);
// 	}


	std::shared_ptr<dasp::KinectGrabber> grabber = std::make_shared<dasp::KinectGrabber>();
	grabber->OpenFile("/home/david/Documents/DataSets/2012-10-12 cogwatch dasp/SlowVelocity/C15_c01_slow.oni");
	grabber->SeekToFrame(100);
	ContinuousSupervoxels sv;
	sv.start(640,480);
	for(int i=0; i<1000; i++) {
		bool ok = grabber->Grab();
		if(!ok) {
			break;
		}
		auto color = grabber->GetLastColor().clone();
		auto depth = grabber->GetLastDepth().clone();
		sv.step(color, depth);
	}

	std::vector<Cluster> clusters = sv.getAllClusters();
	std::cout << "Supervoxel count = " << clusters.size() << std::endl;
	DebugWriteClusters("clusters.tsv", clusters);

	std::cout << "Finished." << std::endl;
//	slimage::gui::WaitForKeypress();

	return 1;
}
