#ifndef INCLUDED_DASP_IO_HPP_
#define INCLUDED_DASP_IO_HPP_

#include "Graph.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace dasp
{
	class Superpixels;

	/** Saves superpixels to a file
	 * Writes rgb (color), xyz (position) and uvw (normal) in that order for each superpixel.
	 * Binary mode: writes raw data (9 x bytes forming a float)
	 * Text mode: writes floats as strings with 3 digits, one superpixel per line
	 */
	void SaveSuperpixels(const Superpixels& superpixels, const std::string& filename, const bool binary=false);

//	void LoadSuperpixels(Superpixels& superpixels, const std::string& filename, bool binary=true);

	void SaveData(const std::vector<std::vector<float>>& data, const std::string& filename, const bool binary=false);

	std::vector<std::vector<float>> LoadData(const std::string& filename, const bool binary=false);

	void SaveGraph(const EdgeWeightGraph& graph, const std::string& filename);
	
	void SaveGraph(const BorderPixelGraph& graph, const std::string& filename);

	struct DaspPoint
	{
		Eigen::Vector2f px;
		Eigen::Vector3f position;
		Eigen::Vector3f color;
		Eigen::Vector3f normal;
	};

	typedef boost::adjacency_list<
		boost::vecS, boost::vecS,
		boost::directedS,
		DaspPoint,
		boost::property<boost::edge_weight_t, float>> DaspGraph;

	DaspGraph LoadDaspGraph(const std::string& fn_dasp, const std::string& fn_graph);


}

#endif
