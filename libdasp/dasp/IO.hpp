#ifndef INCLUDED_DASP_IO_HPP_
#define INCLUDED_DASP_IO_HPP_

#include <string>
#include <vector>

namespace dasp
{
	class Superpixels;

	void SaveSuperpixels(const Superpixels& superpixels, const std::string& filename, const bool binary=true);

//	void LoadSuperpixels(Superpixels& superpixels, const std::string& filename, bool binary=true);

	void SaveData(const std::vector<std::vector<float>>& data, const std::string& filename, const bool binary=true);

	std::vector<std::vector<float>> LoadData(const std::string& filename, const bool binary=true);

}

#endif
