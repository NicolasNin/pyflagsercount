//##############################################################################
//IMPORT LIBRARIES
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <fstream>
#include <vector>
#include <thread>
#include <functional>
#include <deque>
#include <sstream>
#include <cmath>
#include <array>
#include <algorithm>

//##############################################################################
// DEFINITIONS
typedef uint64_t vertex_index_t;
typedef int64_t index_t;
typedef float value_t;

#ifndef PARALLEL_THREADS
#define PARALLEL_THREADS 8
#endif



//##############################################################################
//DIRECTED GRAPH CLASS

bool is_not_in(vertex_index_t el,vertex_index_t l[],unsigned short l_size){
	for (int i=0;i<l_size;i++) {
		if (l[i]==el) return 0;
	}
	return 1;
}
class directed_graphv2_t {
public:
	// The filtration values of the vertices
	vertex_index_t number_of_vertices;
        bool transpose;

	// These are the incidences as a matrix of 64-bit masks
	std::deque<size_t> incidence_outgoing;
	size_t incidence_row_length;

	//we store the transpose for maximal simplices
	//really not necessary but easier this way, we should change this in the future
	std::deque<size_t> incidence_ingoing;
	size_t incidence_col_length;
	// Assume by default that the edge density will be roughly one percent
	directed_graphv2_t(vertex_index_t _number_of_vertices, bool _transpose=false, float density_hint = 0.01)
	    : number_of_vertices(_number_of_vertices), transpose(_transpose), incidence_row_length((_number_of_vertices >> 6) + 1), incidence_col_length((_number_of_vertices >> 6) + 1) {
		incidence_outgoing.resize(incidence_row_length * _number_of_vertices, 0);
		incidence_ingoing.resize(incidence_col_length * _number_of_vertices, 0);
	}

	vertex_index_t vertex_number() const { return number_of_vertices; }

	void add_edge(vertex_index_t v, vertex_index_t w) {
		if(transpose){
                    vertex_index_t v_temp = v;
                    v=w;
                    w=v_temp;
                }
        const size_t ww = w >> 6;
		incidence_outgoing[v * incidence_row_length + ww] |= 1UL << ((w - (ww << 6)));
		
		const size_t vv = v >> 6;
		incidence_ingoing[w * incidence_col_length + vv] |= 1UL << ((v - (vv << 6)));
	}

	bool is_connected_by_an_edge(vertex_index_t from, vertex_index_t to) const {
		const auto t = to >> 6;
		return incidence_outgoing[incidence_row_length * from + t] & (1UL << (to - (t << 6)));
	}

	size_t get_outgoing_chunk(vertex_index_t from, size_t chunk_number) const {
		return incidence_outgoing[incidence_row_length * from + chunk_number];
	}
	size_t get_ingoing_chunk(vertex_index_t from, size_t chunk_number) const {
		return incidence_ingoing[incidence_col_length * from + chunk_number];
	}

	std::vector<vertex_index_t> get_out(vertex_index_t from) const {
		// get all the vertices v coonnected to from as from->v
		std::vector<vertex_index_t> outvertices;
		for (size_t offset = 0; offset < this->incidence_row_length; offset++) {
			size_t bits = this->get_outgoing_chunk(from, offset);

			size_t vertex_offset = offset << 6;
			while (bits > 0) {
				// Get the least significant non-zero bit
				int b = __builtin_ctzl(bits);

				// Unset this bit
				bits &= ~(1UL << b);

				outvertices.push_back(vertex_offset + b);
				}
			}
		return outvertices;
		}

	std::vector<vertex_index_t> get_in(vertex_index_t to) const {
		std::vector<vertex_index_t> invertices;
		for (size_t offset = 0; offset < this->incidence_row_length; offset++) {
			size_t bits = this->get_ingoing_chunk(to, offset);

			size_t vertex_offset = offset << 6;
			while (bits > 0) {
				// Get the least significant non-zero bit
				int b = __builtin_ctzl(bits);

				// Unset this bit
				bits &= ~(1UL << b);

				invertices.push_back(vertex_offset + b);
				}
			}
		return invertices;
		}
	bool vertex_in_between(	vertex_index_t simplex[],unsigned short  simplex_size) const{
		
		//for (int i=0;i<simplex_size;i++) std::cout<< simplex[i]<< " : ";
		//std::cout<<std::endl;
		std::vector<vertex_index_t> candidates;
		std::vector<vertex_index_t> out_start = this->get_out(simplex[0]);
		vertex_index_t end = simplex[simplex_size-1];
		for (int i=0;i<out_start.size();i++){
			if (is_not_in(out_start[i],simplex,simplex_size ) ){
				if (this->is_connected_by_an_edge(out_start[i],end)){ 
					candidates.push_back(out_start[i]);}}
		}

		if (candidates.size()<= 0) return 0;
		// if we found one candidate and its
		if (simplex_size == 2) return 1;

		for (int i=0;i<candidates.size();i++){
			vertex_index_t c = candidates[i];
			bool do_in = 1;
			bool found = 1;
			for (int i=1;i< simplex_size-1;i++){
					if ( do_in && not this->is_connected_by_an_edge(simplex[i],c))
					{
						if (not this->is_connected_by_an_edge(c,simplex[i])) found=0;break;
						do_in=0;
					}
					else
					{
						if (not this->is_connected_by_an_edge(c,simplex[i])) found=0;break;
					}
					
				}
				if (found) return 1;
		}
		// if no candidate is good return false
		return 0;
	}		
};

//##############################################################################
//DIRECTED FLAG COMPLEX CLASS

class directed_flag_complex_max_t {
public:
	const directed_graphv2_t& graph;
	directed_flag_complex_max_t(const directed_graphv2_t& _graph) : graph(_graph) {}

public:
	template <typename Func> void for_each_cell(Func& f, std::vector<vertex_index_t>& do_vertices, std::vector<std::vector<std::vector<vertex_index_t>>>& contain_counts, int min_dimension, int max_dimension = -1) {
		std::array<Func*, 1> fs{&f};
		for_each_cell(fs, do_vertices, contain_counts, min_dimension, max_dimension);
	}

	template <typename Func, size_t number_of_threads>
	void for_each_cell(std::array<Func*, number_of_threads>& fs, std::vector<vertex_index_t>& do_vertices, std::vector<std::vector<std::vector<vertex_index_t>>>& contain_counts, int min_dimension, int max_dimension = -1) {
		if (max_dimension == -1) max_dimension = min_dimension;
		std::thread t[number_of_threads - 1];

		for (size_t index = 0; index < number_of_threads - 1; ++index)
			t[index] = std::thread(&directed_flag_complex_max_t::worker_thread<Func>, this, number_of_threads, index,
			                       fs[index], min_dimension, max_dimension, std::ref(do_vertices), std::ref(contain_counts));

		// Also do work in this thread, namely the last bit
		worker_thread(number_of_threads, number_of_threads - 1, fs[number_of_threads - 1], min_dimension,
		              max_dimension, do_vertices, contain_counts);

		// Wait until all threads stopped
		for (size_t i = 0; i < number_of_threads - 1; ++i) t[i].join();
	}

private:
	template <typename Func>
	void worker_thread(int number_of_threads, int thread_id, Func* f, int min_dimension, int max_dimension,
                       std::vector<vertex_index_t>& do_vertices, std::vector<std::vector<std::vector<vertex_index_t>>>& contain_counts) {
		//const size_t vertices_per_thread = graph.vertex_number() / number_of_threads; //UNUSED ATM

		std::vector<vertex_index_t> first_position_vertices;
		for (size_t index = thread_id; index < do_vertices.size(); index += number_of_threads)
			first_position_vertices.push_back(do_vertices[index]);

		vertex_index_t prefix[max_dimension + 1];

		std::vector<vertex_index_t> possible_prev_vertices;

		do_for_each_cell(f, min_dimension, max_dimension, first_position_vertices, prefix, 0, thread_id, do_vertices.size(), contain_counts, possible_prev_vertices);

		f->done();
	}

	template <typename Func>
	void do_for_each_cell(Func* f, int min_dimension, int max_dimension,
	                      const std::vector<vertex_index_t>& possible_next_vertices, vertex_index_t* prefix,
	                      unsigned short prefix_size, int thread_id, size_t number_of_vertices, std::vector<std::vector<std::vector<vertex_index_t>>>& contain_counts,
						  const std::vector<vertex_index_t>& possible_prev_vertices) {
		// As soon as we have the correct dimension, execute f
		bool is_max = 0;
		
		if (prefix_size>0 && possible_prev_vertices.size()==0 && possible_next_vertices.size()==0 && not graph.vertex_in_between(prefix,prefix_size))	is_max=1;
	
		/*		
		if (is_max == 1){ 
		std::cout<< "IS MAX ?"<<std::endl;
		std::cout<< "prefix size "<<prefix_size  <<std::endl;
		
			for (auto i =0;i<prefix_size;i++){
				std::cout<< prefix[i]<< " " ; 
				}
				std::cout<<std::endl;
			}
		*/
		if (prefix_size >= min_dimension + 1) {
			 (*f)(prefix, prefix_size,is_max); 
			 }
		

		// If this is the last dimension we are interested in, exit this branch
		if (prefix_size == max_dimension + 1) return;


        for (auto vertex : possible_next_vertices) {
			// We can write the cell given by taking the current vertex as the maximal element
			prefix[prefix_size] = vertex;

			// And compute the next elements
			std::vector<vertex_index_t> new_possible_vertices;
			std::vector<vertex_index_t> new_possible_prev_vertices;

			if (prefix_size > 0) {
				for (auto v : possible_next_vertices) {
					if (vertex != v && graph.is_connected_by_an_edge(vertex, v)) new_possible_vertices.push_back(v);

					for (auto vprev : possible_prev_vertices) 
					if (vertex != vprev && graph.is_connected_by_an_edge(vprev, vertex)) new_possible_prev_vertices.push_back(vprev);

				}
			} else {
				// Get outgoing vertices of v in chunks of 64
				for (size_t offset = 0; offset < graph.incidence_row_length; offset++) {
					size_t bits = graph.get_outgoing_chunk(vertex, offset);
					// we repeat everything for ingoing
					size_t inbits = graph.get_ingoing_chunk(vertex, offset);

					size_t vertex_offset = offset << 6;
					while (bits > 0) {
						// Get the least significant non-zero bit
						int b = __builtin_ctzl(bits);
						// Unset this bit
						bits &= ~(1UL << b);
						new_possible_vertices.push_back(vertex_offset + b);
					}
					while (inbits > 0) {
						// Get the least significant non-zero bit
						int b = __builtin_ctzl(inbits);
						// Unset this bit
						inbits &= ~(1UL << b);
						new_possible_prev_vertices.push_back(vertex_offset + b);
					}
				}
			}

            do_for_each_cell(f, min_dimension, max_dimension, new_possible_vertices, prefix, prefix_size + 1, thread_id, number_of_vertices, contain_counts, new_possible_prev_vertices);
		}
	}
};


//##############################################################################
//CELL COUNTER STRUCT

// we modify a tiny bit the  cell_counter_t previous struct so that there is two counter one for maximal and one for the rest
struct cell_counter_max_t {
	void done() {}
	void operator()(vertex_index_t* first_vertex, int size, bool is_max) {
		// Add (-1)^size to the Euler characteristic
		if (size & 1)
			ec++;
		else
			ec--;

		if (cell_counts.size() < size) { cell_counts.resize(size, 0); }
		cell_counts[size - 1]++;

		if (is_max){
		if (cell_max_counts.size() < size) { cell_max_counts.resize(size, 0); }
		cell_max_counts[size - 1]++;
		}
	}

	int64_t euler_characteristic() const { return ec; }
	std::vector<size_t> cell_count() const { return cell_counts; }
	std::vector<size_t> cell_max_count() const { return cell_max_counts; }

private:
	int64_t ec = 0;
	std::vector<size_t> cell_counts;
	std::vector<size_t> cell_max_counts;

};
//##############################################################################
//COUNT CELL FUNCTION

std::vector<std::vector<size_t>> count_cells_max(directed_graphv2_t& graph) {
	directed_flag_complex_max_t complex(graph);

   std::vector<vertex_index_t> do_vertices;
   for(int i = 0; i < graph.vertex_number(); i++){ do_vertices.push_back(i); }

	//TODO REMOVE CONTAINS everywhere
    std::vector<std::vector<std::vector<vertex_index_t>>> contain_counts(PARALLEL_THREADS,
			                                                     std::vector<std::vector<vertex_index_t>>(graph.vertex_number(),
																														std::vector<vertex_index_t>(0)));

	std::array<cell_counter_max_t*, PARALLEL_THREADS> cell_counter;
	for (int i = 0; i < PARALLEL_THREADS; i++)
		cell_counter[i] = new cell_counter_max_t();

		complex.for_each_cell(cell_counter, do_vertices, contain_counts, 0, 10000);
	
	std::vector<size_t> total_counts;
	std::vector<size_t> total_max_counts;
	for (int i=0;i<PARALLEL_THREADS;i++){
		std::vector<size_t> counti = cell_counter[i]->cell_count();
		if (total_counts.size()<counti.size()) { total_counts.resize(counti.size(), 0); }
		for (int j=0;j<counti.size();j++){ total_counts[j]+=counti[j]; }

		std::vector<size_t> count_maxi = cell_counter[i]->cell_max_count();
		if (total_max_counts.size()<count_maxi.size()) { total_max_counts.resize(count_maxi.size(), 0); }
		for (int j=0;j<count_maxi.size();j++){ total_max_counts[j]+=count_maxi[j]; }

	}

	std::cout<<"* Total SIMPLICES"<<std::endl;
	for (int i=0; i<PARALLEL_THREADS; i++){
		std::cout<< i<< ": ";
		std::vector<size_t> counti = cell_counter[i]->cell_count();
		for (int j=0;j<counti.size();j++)
			std::cout<< counti[j]<< " ";
		std::cout<< std::endl;
	}
	for (int i=0;i<total_counts.size();i++){std::cout<< total_counts[i]<<" ";}

	std::cout<<std::endl<<"* MAXIMAL SIMPLICES"<<std::endl;
		for (int i=0; i<PARALLEL_THREADS; i++){
		std::cout<< i<< ": ";
		std::vector<size_t> counti = cell_counter[i]->cell_max_count();
		for (int j=0;j<counti.size();j++)
			std::cout<< counti[j]<< " ";
		std::cout<< std::endl;
	}
	for (int i=0;i<total_max_counts.size();i++){std::cout<< total_max_counts[i]<<" ";}
	std::cout<<std::endl;
	std::vector<std::vector<size_t>> ret;
	ret.push_back(total_counts);
	ret.push_back(total_max_counts);
	
	return ret;
}
/*

int main(int argc, char** argv){


	int N=3;
	auto graph = directed_graphv2_t(N);

	//std::vector<vertex_index_t> row{0,0,1};
	//std::vector<vertex_index_t> col{1,2,2};
	// undirected triangle
	std::vector<vertex_index_t> row{0,1,2};
	std::vector<vertex_index_t> col{1,2,0};
	std::vector<vertex_index_t> data{1,8,1};
	for (int i=0;i<3;i++)
		graph.add_edge(row[i],col[i]);

	std::cout<< "Matrix"<<std::endl;
	for (vertex_index_t i=0;i<3;i++){
		for (vertex_index_t j=0;j<3;j++){
			if (graph.is_connected_by_an_edge(i,j))
			{std::cout<< 1<< " ";}
			else{std::cout<< "0 ";}
			}
	std::cout<<std::endl;
	}
	std::cout<< "------------------"<<std::endl;
	std::vector<std::vector<size_t>> a =  count_cells_max(graph);

	
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
	}
}
*/