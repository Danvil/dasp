/*
 * Segmentation.cpp
 *
 *  Created on: Mar 26, 2012
 *      Author: david
 */

//#define SEGS_DBG_SHOWGUI
//#define SEGS_DBG_CREATE_EV_IMAGES
//#define SEGS_DBG_PRINT
//#define SEGS_VERBOSE

namespace dasp
{

#ifdef SEGS_DBG_SHOWGUI
	{
		slimage::Image3ub vis = clusters.color_raw.clone();
		dasp::plots::PlotEdges(vis, clusters.ComputeLabels(), slimage::Pixel3ub{{255,255,255}},1);
		slimage::gui::Show("color", vis);
	}
	{
		slimage::Image3ub vis = dasp::plots::PlotClusters(clusters, dasp::plots::ClusterPoints, dasp::plots::Color);
		dasp::plots::PlotEdges(vis, clusters.ComputeLabels(), slimage::Pixel3ub{{255,255,255}},1);
		slimage::gui::Show("clusters", vis);
	}

#endif

#ifdef SEGS_VERBOSE
	std::cout << "SpectralSegmentation: n = " << n << std::endl;
#endif

#ifdef SEGS_VERBOSE
	std::cout << "Edve connectivity: min=" << edge_connectivity.minCoeff() << ", max=" << edge_connectivity.maxCoeff() << std::endl;
#endif

#ifdef SEGS_DBG_SHOWGUI
	{
		slimage::Image3ub edge_connectivity_img(clusters.width(), clusters.height(), slimage::Pixel3ub{{0,0,0}});
		for(unsigned int eid=0; eid<Gnb.edges.size(); eid++) {
			for(unsigned int pid : border_pixels[eid]) {
				edge_connectivity_img[pid] = dasp::plots::IntensityColor(static_cast<float>(edge_connectivity[eid]), 0.0f, 1.0f);
			}
		}
		slimage::gui::Show("edge_connectivity", edge_connectivity_img);
	}
#endif

#ifdef SEGS_VERBOSE
	std::cout << "Eigenvalues = " << result_ew.topRows(n_used_ew + 1).transpose() << std::endl;
#endif
#ifdef SEGS_DBG_CREATE_EV_IMAGES
	{
		cSegmentationDebug.clear();
		// create image from eigenvectors (omit first)
		for(unsigned int k=0; k<std::min(cNEV,3u); k++) {
			// get k-th eigenvector
			Vec ev = result_ev.col(k + 1);
			// convert to plotable values
			std::vector<unsigned char> ev_ub(n);
			for(unsigned int i=0; i<n; i++) {
				float v = 0.5f + 2.0f*ev[i];
				ev_ub[i] = static_cast<unsigned char>(std::min(255, std::max(0, static_cast<int>(255.0f * v))));
			}
			// write to image
			slimage::Image3ub img(clusters.width(), clusters.height(), slimage::Pixel3ub{{255,0,0}});
			clusters.ForPixelClusters([&img,&ev_ub](unsigned int cid, const dasp::Cluster& c, unsigned int pid, const dasp::Point& p) {
				unsigned char v = ev_ub[cid];
				img[pid] = slimage::Pixel3ub{{v,v,v}};
			});
			cSegmentationDebug.push_back(img);
#ifdef SEGS_DBG_SHOWGUI
			slimage::gui::Show((boost::format("ev_%02d") % (k+1)).str(), img);
#endif
		}
	}	// DEBUG
#endif

#ifdef SEGS_VERBOSE
		std::cout << "w=" << w << " e_k.maxCoeff()=" << e_k.maxCoeff() << std::endl;
#endif

#ifdef SEGS_DBG_PRINT
		{
			std::ofstream ofs((boost::format("/tmp/edge_weights_%03d.txt") % k).str());
			for(unsigned int i=0; i<e_k.rows(); i++) {
				ofs << e_k[i] << std::endl;
			}
		}
#endif

#ifdef SEGS_DBG_PRINT
	{
		std::ofstream ofs("/tmp/edge_weights_sum.txt");
		for(unsigned int i=0; i<edge_weight.rows(); i++) {
			ofs << edge_weight[i] << std::endl;
		}
	}
#endif


#ifdef SEGS_VERBOSE
	std::cout << "Edge weights: min=" << edge_weight.minCoeff() << ", max=" << edge_weight.maxCoeff() << std::endl;
	std::cout << "Edge weights: median=" << edge_weight[edge_weight.rows()/2] << std::endl;
#endif

#ifdef SEGS_DBG_SHOWGUI
	slimage::Image1ub boundaries = CreateBorderImage(clusters.width(), clusters.height(), graph); // FIXME <- fuse border_pixels
	slimage::gui::Show("boundaries", segs.boundaries, 0.03f);
#endif

#ifdef SEGS_DBG_SHOWGUI
	slimage::gui::WaitForKeypress();
#endif

}
