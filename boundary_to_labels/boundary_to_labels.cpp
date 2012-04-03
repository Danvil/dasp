#include <Slimage/IO.hpp>
#include <Slimage/Slimage.hpp>
#include <Slimage/Convert.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <iostream>
#include <string>

slimage::Image1ui16 BoundaryToLabels(const slimage::Image1ub& img_bnds)
{
	slimage::Image1i temp(img_bnds.dimensions());

	int label = 1;

	// do 1d flood fill for each line
	for(unsigned int y=0; y<img_bnds.height(); y++) {
		bool last_was_border = true;
		for(unsigned int x=0; x<img_bnds.width(); x++) {
			bool is_border = img_bnds(x,y) > 0;
			if(is_border) {
				if(!last_was_border)
					label++;
				temp(x,y) = 0;
			}
			else {
				temp(x,y) = label;
			}
			last_was_border = is_border;
		}
	}

	// now combine neighbouring lines
	std::map<int, std::set<int> > neighbours;
	for(unsigned int y=0; y+1<img_bnds.height(); y++) {
		for(unsigned int x=0; x<img_bnds.width(); x++) {
			int label = temp(x, y);
			int label_down = temp(x, y+1);
			if(label != 0 && label_down != 0) {
				neighbours[label].insert(label_down);
				neighbours[label_down].insert(label);
			}
		}
	}

	// find connected components
	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS> MyGraph;
	MyGraph graph;
	std::map<int, MyGraph::vertex_descriptor> label_to_vid;
	for(auto p : neighbours) {
		MyGraph::vertex_descriptor vid = boost::add_vertex(graph);
		label_to_vid[p.first] = vid;
	}
	for(auto p : neighbours) {
		for(int v : p.second) {
			boost::add_edge(label_to_vid[p.first], label_to_vid[v], graph);
		}
	}
	// compute connected components
	typedef std::map<MyGraph::vertex_descriptor, MyGraph::vertices_size_type> component_type;
	component_type component;
	boost::associative_property_map< component_type > component_map(component);
	unsigned int segment_count = boost::connected_components(graph, component_map);
	std::cout << "Number of components: " << segment_count << std::endl;

	// write results
	slimage::Image1ui16 result(img_bnds.dimensions());
	for(unsigned int y=0; y<img_bnds.height(); y++) {
		for(unsigned int x=0; x<img_bnds.width(); x++) {
			int q = temp(x,y);
			if(q == 0) {
				result(x,y) = static_cast<uint16_t>(0);
			}
			else {
				result(x,y) = static_cast<uint16_t>(component[label_to_vid[q]] + 1);
			}
		}
	}

	while(true) {
		unsigned int cnt_still_border = 0;
		for(unsigned int y=0; y<img_bnds.height(); y++) {
			for(unsigned int x=0; x<img_bnds.width(); x++) {
				if(result(x,y) == 0) {
					std::set<int> nbs;
					if(1 <= x) {
						nbs.insert(result(x-1,y));
					}
					if(x+1 < img_bnds.width()) {
						nbs.insert(result(x+1,y));
					}
					if(1 <= y) {
						nbs.insert(result(x,y-1));
					}
					if(y+1 < img_bnds.width()) {
						nbs.insert(result(x,y+1));
					}
					nbs.erase(0);
					if(nbs.size() == 0) {
						cnt_still_border++;
					}
					else {
						// take some
						result(x,y) = *(nbs.begin());
					}
				}
			}
		}
		if(cnt_still_border == 0) {
			break;
		}
	}

	return result;
}

int main(int argc, char** argv)
{
	std::string fn_in = argv[1];
	std::string fn_out = argv[2];

	slimage::Image1ub img_bnds = slimage::Pick<unsigned char>(slimage::Load(fn_in), 0);

	slimage::Image1ui16 img_labels = BoundaryToLabels(img_bnds);
	slimage::Save(img_labels, fn_out);

	return 0;
}
