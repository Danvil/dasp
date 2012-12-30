#include <dasv.hpp>
#include <graphseg/Rendering.hpp>
#include <Candy/System/GlutSystem.h>
#include <boost/program_options.hpp>

int main(int argc, char** argv)
{
	std::string fn_clusters = "/tmp/dasv/00020_clusters.tsv";
	std::string fn_edges = "/tmp/dasv/00020_edges.tsv";

	namespace po = boost::program_options;
	po::options_description desc;
	desc.add_options()
		("help", "produce help message")
		("vertices", po::value<std::string>(&fn_clusters)->required(), "filename to vertices file")
		("edges", po::value<std::string>(&fn_edges)->required(), "filename to edges file")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
	if(vm.count("help")) {
		std::cerr << desc << std::endl;
		return 1;
	}

	dasv::ClusterGraph graph = dasv::IOReadGraph(fn_clusters, fn_edges);

	Candy::ViewPtr view = Candy::View::FactorDefaultPerspectiveView();
	Candy::ScenePtr scene = view->getScene();

	boost::shared_ptr<Candy::IRenderable> renderling(new Candy::ObjectRenderling(
		[&graph]() {
			graphseg::RenderEdges3DVcol(graph,
				[&graph](const dasv::ClusterGraph::vertex_descriptor& vid) { // returns vertex coordinate
					const dasv::Cluster& c = graph[vid];
					return Eigen::Vector3f{ 0.01f*c.pixel.x(), 0.01f*c.pixel.y(), 0.03f*static_cast<float>(c.time) };
					// return c.position;
				},
				[&graph](const dasv::ClusterGraph::vertex_descriptor& vid) { // returns vertex color
					return graph[vid].color;
				}
			);
		}
	));
	scene->addItem(renderling);

	Candy::EnginePtr engine(new Candy::Engine(view));
	engine->setClearColor(Danvil::Color::DarkGrey);
	scene->setShowCoordinateCross(true);

	Candy::GlutSystem sys(engine, "dasp/dasv graph visualization");
	return sys.main(argc, argv);
}
