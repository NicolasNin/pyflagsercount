#include "flagser-count-maximal.cpp"
int main(int argc, char** argv){



	//std::vector<vertex_index_t> row{0,0,1};
	//std::vector<vertex_index_t> col{1,2,2};
	// undirected triangle
//	std::vector<vertex_index_t> row{0,1,2};
//	std::vector<vertex_index_t> col{1,2,0};
   // std::vector<vertex_index_t> row{0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 4, 4};
   // std::vector<vertex_index_t> col{1, 2, 3, 4, 0, 3, 4, 0, 1, 0, 1, 2};
   // std::vector<vertex_index_t> row{0, 0, 1, 1, 1, 2, 2, 2, 3};
    //std::vector<vertex_index_t> col{1, 2, 0, 2, 3, 0, 1, 3, 0};

  std::vector<vertex_index_t> row   {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4};
 std::vector<vertex_index_t> col{1, 3, 4, 2, 3, 4, 0, 1, 3, 4, 2, 4, 0, 2, 3};

    int N=5;
	auto graph = directed_graphv2_t(N);

	for (int i=0;i<row.size();i++)
		graph.add_edge(row[i],col[i]);

	std::cout<< "Matrix"<<std::endl;
	for (vertex_index_t i=0;i<N;i++){
		for (vertex_index_t j=0;j<N;j++){
			if (graph.is_connected_by_an_edge(i,j))
			{std::cout<< 1<< " ";}
			else{std::cout<< "0 ";}
			}
	std::cout<<std::endl;
	}
	std::cout<< "------------------"<<std::endl;
	std::vector<std::vector<size_t>> a =  count_cells_max(graph);
    vertex_index_t t[4]={2,4,0,3};
    
	std::cout<<" in between "<< std::endl<<graph.vertex_in_between(t,4)<<std::endl;
	/*
	auto out0=graph.get_in(2);
	std::cout<<"OUT 0" <<std::endl;
	for (int i=0;i<out0.size();i++)
		std::cout<< out0[i]<< " ";
	std::cout<< std::endl;
	
	directed_flag_complex_max_t complex(graph);
	std::vector<vertex_index_t> do_vertices;
	//TODO REMOVE CONTAINS everywhere
   for(int i = 0; i < graph.vertex_number(); i++){ do_vertices.push_back(i); }
	std::vector<std::vector<std::vector<vertex_index_t>>> contain_counts(PARALLEL_THREADS,
			                                                     std::vector<std::vector<vertex_index_t>>(graph.vertex_number(),
																														std::vector<vertex_index_t>(0)));

	std::array<cell_counter_max_t*, PARALLEL_THREADS> cell_counter;
	for (int i = 0; i < PARALLEL_THREADS; i++)
		cell_counter[i] = new cell_counter_max_t();

	
	complex.for_each_cell(cell_counter, do_vertices, contain_counts, 0, 10000);
	
	std::cout<<"* Total SIMPLICES"<<std::endl;
	for (int i=0; i<PARALLEL_THREADS; i++){
		std::cout<< i<< ": ";
		std::vector<size_t> counti = cell_counter[i]->cell_count();
		for (int j=0;j<counti.size();j++)
			std::cout<< counti[j]<< " ";
		std::cout<< std::endl;
	}
	std::cout<<"* MAXIMAL SIMPLICES"<<std::endl;
		for (int i=0; i<PARALLEL_THREADS; i++){
		std::cout<< i<< ": ";
		std::vector<size_t> counti = cell_counter[i]->cell_max_count();
		for (int j=0;j<counti.size();j++)
			std::cout<< counti[j]<< " ";
		std::cout<< std::endl;
	}*/
}