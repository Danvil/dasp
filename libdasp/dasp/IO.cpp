#include "IO.hpp"
#include "Superpixels.hpp"
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace dasp
{
	void SaveSuperpixels(const Superpixels& superpixels, const std::string& filename, const bool binary)
	{
		std::vector<std::vector<float>> data;
		data.reserve(superpixels.cluster.size());
		for(const Cluster& c : superpixels.cluster) {
			const Point& p = c.center;
			Eigen::Vector3f n = p.computeNormal();
			std::vector<float> v {
				p.world.x(), p.world.y(), p.world.z(),
				p.color.x(), p.color.y(), p.color.z(),
				n.x(), n.y(), n.z()
			};
			data.push_back(v);
		}
		SaveData(data, filename, binary);
	}

	// void LoadSuperpixels(Superpixels& superpixels, const std::string& filename, bool binary)
	// {
	// 	throw "Not implemented";
	// 	std::vector<std::vector<float>> data = LoadData(filename, binary);
	// 	superpixels.cluster.clear();
	// 	superpixels.cluster.reserve(data.size());
	// 	for(const std::vector<float>& v : data) {
	// 		//
	// 	}
	// }

	void SaveData(const std::vector<std::vector<float>>& data, const std::string& filename, const bool binary)
	{
		std::ofstream ofs;
		if(binary) {
			ofs.open(filename, std::fstream::out | std::fstream::binary);
		}
		else {
			ofs.open(filename, std::fstream::out);
		}
		if(!ofs.is_open()) {
			// FIXME throw proper exception
			throw "Could not open file";
		}
		if(data.size() == 0) {
			return;
		}
		ofs.precision(4); // 4 is sub mm for position - 3 would be enough for color and normal
		const unsigned int size = data[0].size();
		for(const std::vector<float>& v : data) {
			// if(v.size() != size) {
			// 	throw "Invalid data!";
			// }
			const float* p = &v[0];
			if(binary) {
				ofs.write(reinterpret_cast<const char*>(p), size*sizeof(float));
			}
			else {
				for(unsigned int i=0; i<size; i++) {
					ofs << p[i];
					if(i == size-1)
						ofs << "\n";
					else
						ofs << "\t";
				}
			}
		}
	}

	std::vector<std::vector<float>> LoadData(const std::string& filename, const bool binary)
	{
		std::ifstream ifs;
		if(binary) {
			ifs.open(filename, std::fstream::in | std::fstream::binary);
		}
		else {
			ifs.open(filename, std::fstream::in);
		}
		if(!ifs.is_open()) {
			// FIXME throw exception
			return {};
		}
		boost::escaped_list_separator<char> separator('\\','\t','\"');
		typedef boost::tokenizer<decltype(separator)> Tokenizer;
		std::vector<std::vector<float>> data;
		std::string line;
		std::vector<float> v;
		while(std::getline(ifs,line))
		{
			v.clear();
			Tokenizer tok(line, separator);
			for(Tokenizer::iterator beg=tok.begin(); beg!=tok.end(); ++beg) {
				float x = boost::lexical_cast<float>(*beg);
				v.push_back(x);
			}
			data.push_back(v);
		}
		return data;
	}

}

