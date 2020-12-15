#include "filtered_flagser-count.cpp"

int main(int argc, char** argv) {
	std::cout<<"test"<<std::endl;
    auto graph = directed_graph_t(10);

		std::vector<vertex_index_t> row{0,0,1};
		std::vector<vertex_index_t> col{1,2,2};
		std::vector<vertex_index_t> data{1,8,1};


		auto graph2 = compressed_directed_graph_t(3,false);

		for (int i=0;i<3;i++)
			graph2.add_edge_weighted(row[i],col[i],data[i]);

		// graph checks
		//###############
		std::cout<<"number of vertex "<<graph2.vertex_number()<<std::endl;
		std::cout<<"G[0 1] "<<graph2.is_connected_by_an_edge(0,1)<<std::endl;
		std::cout<<"Max filtration value "<<graph2.max_filtration()<<std::endl;
		std::vector<vertex_index_t> new_possible_vertices;
		hash_map* out_neigh = graph2.get_outgoing_chunk(0);
	 for(auto iter = out_neigh->begin(); iter != out_neigh->end(); ++iter){
		 new_possible_vertices.push_back(iter->first);}
		for (int i=0;i<new_possible_vertices.size();i++)
			std::cout<< "  " <<new_possible_vertices[i];
		std::cout<<std::endl;

		for (vertex_index_t i=0;i<3;i++){
			for (vertex_index_t j=0;j<3;j++){
				if (graph2.is_connected_by_an_edge(i,j))
				{std::cout<< graph2.edge_value(i,j)<< " ";}else{std::cout<< "0 ";}
			}
			std::cout<<std::endl;
		}
		// END GRAPH CHECKS
			auto filtration_counts = filtered_count_cells(graph2);
			for (int i=0;i<filtration_counts.size();i++){
				std::cout<<i<< ":";
				for(int j=0;j<filtration_counts[i].size();j++)
					std::cout<< filtration_counts[i][j]<< " ";
				std::cout<<std::endl;
			}
		}
