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
#include <google/dense_hash_map>

//##############################################################################
// DEFINITIONS
typedef uint64_t vertex_index_t;
typedef int64_t index_t;
typedef float value_t;

#define PARALLEL_THREADS 72

class hash_map : public google::dense_hash_map<vertex_index_t, vertex_index_t> {
  public:
	  inline void reserve(size_t hint) { this->resize(hint); }
};

//##############################################################################
//DIRECTED GRAPH CLASS

class directed_graph_t {
public:
	// The filtration values of the vertices
	vertex_index_t number_of_vertices;
        bool transpose;

	// These are the incidences as a matrix of 64-bit masks
	std::deque<size_t> incidence_outgoing;
	size_t incidence_row_length;

	// Assume by default that the edge density will be roughly one percent
	directed_graph_t(vertex_index_t _number_of_vertices, bool _transpose=false, float density_hint = 0.01)
	    : number_of_vertices(_number_of_vertices), transpose(_transpose), incidence_row_length((_number_of_vertices >> 6) + 1) {
		incidence_outgoing.resize(incidence_row_length * _number_of_vertices, 0);
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
	}

	bool is_connected_by_an_edge(vertex_index_t from, vertex_index_t to) const {
		const auto t = to >> 6;
		return incidence_outgoing[incidence_row_length * from + t] & (1UL << (to - (t << 6)));
	}

	size_t get_outgoing_chunk(vertex_index_t from, size_t chunk_number) const {
		return incidence_outgoing[incidence_row_length * from + chunk_number];
	}
};
//#############################################################################
//Compressed Class
//Stores the adjacency matrix as a vector of dense_hash_maps, one for each vertex
class compressed_directed_graph_t : public directed_graph_t {
public:
	std::vector<hash_map> incidence_outgoing;
	compressed_directed_graph_t(vertex_index_t _number_of_vertices, bool _transpose=false)
	    : directed_graph_t{ _transpose } {
    set_number_of_vertices(_number_of_vertices);
	}
	vertex_index_t max_filtration_value = 0;
	vertex_index_t max_filtration()  { return max_filtration_value; }

	virtual void set_number_of_vertices(vertex_index_t _number_of_vertices){
		if(number_of_vertices != _number_of_vertices){
	    number_of_vertices = _number_of_vertices;
			incidence_outgoing.clear();
			for(int i=0; i < _number_of_vertices; i++){
					incidence_outgoing.push_back(hash_map());
					incidence_outgoing[i].set_empty_key(std::numeric_limits<vertex_index_t>::max());
			}
	  }
	}

	virtual void add_edge(vertex_index_t v, vertex_index_t w) {
		if(v >= number_of_vertices || w >= number_of_vertices){
			std::cerr << "ERROR: Edge " << v << " " << w << " can't exist, as largest vertex id is " << number_of_vertices-1 << std::endl;
			exit(-1);
		}
	  if(transpose){ incidence_outgoing[w][v] = 1; }
    else{ incidence_outgoing[v][w] = 1; }
  }
	virtual void add_edge_weighted(vertex_index_t v, vertex_index_t w,vertex_index_t edge_value) {
		max_filtration_value=std::max(max_filtration_value,edge_value);
		if(v >= number_of_vertices || w >= number_of_vertices){
			std::cerr << "ERROR: Edge " << v << " " << w << " can't exist, as largest vertex id is " << number_of_vertices-1 << std::endl;
			exit(-1);
		}
	  if(transpose){ incidence_outgoing[w][v] = edge_value; }
    else{ incidence_outgoing[v][w] = edge_value; }
  }
	virtual bool is_connected_by_an_edge(vertex_index_t from, vertex_index_t to) const{
    return (incidence_outgoing[from].find(to) != incidence_outgoing[from].end());
	}
	virtual vertex_index_t edge_value(vertex_index_t from, vertex_index_t to) {
    return incidence_outgoing[from][to] ;
	}
	//returns out neighbours
	virtual hash_map*  get_outgoing_chunk(vertex_index_t from) {
		return  &incidence_outgoing[from];
	}
};
//##############################################################################
//DIRECTED FLAG COMPLEX CLASS

class directed_flag_complex_t {
public:
	const directed_graph_t& graph;
	directed_flag_complex_t(const directed_graph_t& _graph) : graph(_graph) {}

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
			t[index] = std::thread(&directed_flag_complex_t::worker_thread<Func>, this, number_of_threads, index,
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
		const size_t vertices_per_thread = graph.vertex_number() / number_of_threads;

		std::vector<vertex_index_t> first_position_vertices;
		for (size_t index = thread_id; index < do_vertices.size(); index += number_of_threads)
			first_position_vertices.push_back(do_vertices[index]);

		vertex_index_t prefix[max_dimension + 1];

		do_for_each_cell(f, min_dimension, max_dimension, first_position_vertices, prefix, 0, thread_id, do_vertices.size(), contain_counts);

		f->done();
	}

	template <typename Func>
	void do_for_each_cell(Func* f, int min_dimension, int max_dimension,
	                      const std::vector<vertex_index_t>& possible_next_vertices, vertex_index_t* prefix,
	                      unsigned short prefix_size, int thread_id, size_t number_of_vertices, std::vector<std::vector<std::vector<vertex_index_t>>>& contain_counts) {
		// As soon as we have the correct dimension, execute f
		if (prefix_size >= min_dimension + 1) { (*f)(prefix, prefix_size); }
        for(int i = 0; i < prefix_size; i++){
            while(contain_counts[thread_id][prefix[i]].size() < prefix_size){
                contain_counts[thread_id][prefix[i]].push_back(0);
            }
            contain_counts[thread_id][prefix[i]][prefix_size-1]++;
        }

		// If this is the last dimension we are interested in, exit this branch
		if (prefix_size == max_dimension + 1) return;

        for (auto vertex : possible_next_vertices) {
			// We can write the cell given by taking the current vertex as the maximal element
			prefix[prefix_size] = vertex;

			// And compute the next elements
			std::vector<vertex_index_t> new_possible_vertices;
			if (prefix_size > 0) {
				for (auto v : possible_next_vertices) {
					if (vertex != v && graph.is_connected_by_an_edge(vertex, v)) new_possible_vertices.push_back(v);
				}
			} else {
				get_new_possible_vertex(vertex, new_possible_vertices);}

            do_for_each_cell(f, min_dimension, max_dimension, new_possible_vertices, prefix, prefix_size + 1, thread_id, number_of_vertices, contain_counts);
		}
	}
	virtual void get_new_possible_vertex(vertex_index_t vertex, std::vector<vertex_index_t>& new_possible_vertices) {
		for (size_t offset = 0; offset < graph.incidence_row_length; offset++) {
			size_t bits = graph.get_outgoing_chunk(vertex, offset);

			size_t vertex_offset = offset << 6;
			while (bits > 0) {
				int b = __builtin_ctzl(bits);  // Get the least significant non-zero bit
				bits &= ~(1UL << b);           // Unset this bit
				new_possible_vertices.push_back(vertex_offset + b);
			}
		}
	}

};
//##############################################################################
//FILTERED DIRECTED FLAG COMPLEX CLASS

class filtered_directed_flag_complex_t {
public:
	 compressed_directed_graph_t& graph;
	filtered_directed_flag_complex_t( compressed_directed_graph_t& _graph) : graph(_graph) {}
	/*const compressed_directed_graph_t& graph;
	filtered_directed_flag_complex_t(const compressed_directed_graph_t& _graph) : graph(_graph) {}
*/
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
			t[index] = std::thread(&filtered_directed_flag_complex_t::worker_thread<Func>, this, number_of_threads, index,
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
		const size_t vertices_per_thread = graph.vertex_number() / number_of_threads;

		std::vector<vertex_index_t> first_position_vertices;
		for (size_t index = thread_id; index < do_vertices.size(); index += number_of_threads)
			first_position_vertices.push_back(do_vertices[index]);

		vertex_index_t prefix[max_dimension + 1];

		do_for_each_cell(f, min_dimension, max_dimension, first_position_vertices, prefix, 0, thread_id, do_vertices.size(),contain_counts, 0);

		f->done();
	}

	template <typename Func>
	void do_for_each_cell(Func* f, int min_dimension, int max_dimension,
	                      const std::vector<vertex_index_t>& possible_next_vertices, vertex_index_t* prefix,
	                      unsigned short prefix_size, int thread_id, size_t number_of_vertices, std::vector<std::vector<std::vector<vertex_index_t>>>& contain_counts,
											vertex_index_t filtration_value) {
// I WILL OVERRIDE THE CONTAINS COUNTS TO STORE THE FILTRATION COUNTS
		// As soon as we have the correct dimension, execute f
		if (prefix_size >= min_dimension + 1) { (*f)(prefix, prefix_size); }

		if (prefix_size>=1){
            while(contain_counts[thread_id][filtration_value].size() < prefix_size){
                contain_counts[thread_id][filtration_value].push_back(0);
            }
            contain_counts[thread_id][filtration_value][prefix_size-1]++;
		}

		// If this is the last dimension we are interested in, exit this branch
		if (prefix_size == max_dimension + 1) return;
        for (auto vertex : possible_next_vertices) {
			// We can write the cell given by taking the current vertex as the maximal element
			prefix[prefix_size] = vertex;
			// update the new filtration_value
			vertex_index_t new_filtration_value=filtration_value;

			for (auto i=0;i<prefix_size;i++)
				new_filtration_value=std::max(new_filtration_value, graph.edge_value(prefix[i],vertex));
			// And compute the next elements
			std::vector<vertex_index_t> new_possible_vertices;
			if (prefix_size > 0) {
				for (auto v : possible_next_vertices) {
					if (vertex != v && graph.is_connected_by_an_edge(vertex, v)) new_possible_vertices.push_back(v);
				}
			} else {
				get_new_possible_vertex(vertex, new_possible_vertices);}

            do_for_each_cell(f, min_dimension, max_dimension, new_possible_vertices, prefix, prefix_size + 1, thread_id, number_of_vertices, contain_counts,new_filtration_value);
		}
	}

	virtual void get_new_possible_vertex(vertex_index_t vertex, std::vector<vertex_index_t>& new_possible_vertices) {
		 hash_map* out_neigh = graph.get_outgoing_chunk(vertex);
		for(auto iter = out_neigh->begin(); iter != out_neigh->end(); ++iter){
			new_possible_vertices.push_back(iter->first);
		}
	}

};

//##############################################################################
//CELL COUNTER STRUCT

struct cell_counter_t {
	void done() {}
	void operator()(vertex_index_t* first_vertex, int size) {
		// Add (-1)^size to the Euler characteristic
		if (size & 1)
			ec++;
		else
			ec--;

		if (cell_counts.size() < size) { cell_counts.resize(size, 0); }
		cell_counts[size - 1]++;
	}

	int64_t euler_characteristic() const { return ec; }
	std::vector<size_t> cell_count() const { return cell_counts; }

private:
	int64_t ec = 0;
	std::vector<size_t> cell_counts;
};

//##############################################################################
//COUNT CELL FUNCTION

std::vector<std::vector<vertex_index_t>> count_cells(directed_graph_t& graph) {
	directed_flag_complex_t complex(graph);


   std::vector<vertex_index_t> do_vertices;
   for(int i = 0; i < graph.vertex_number(); i++){ do_vertices.push_back(i); }


    std::vector<std::vector<std::vector<vertex_index_t>>> contain_counts(PARALLEL_THREADS,
			                                                     std::vector<std::vector<vertex_index_t>>(graph.vertex_number(),
																														std::vector<vertex_index_t>(0)));

	std::array<cell_counter_t*, PARALLEL_THREADS> cell_counter;
	for (int i = 0; i < PARALLEL_THREADS; i++)
		cell_counter[i] = new cell_counter_t();

		complex.for_each_cell(cell_counter, do_vertices, contain_counts, 0, 10000);


    for(int i = 1; i < contain_counts.size(); i++){
        for(int j = 0; j < contain_counts[i].size(); j++){
            while(contain_counts[0][j].size() < contain_counts[i][j].size()){
                contain_counts[0][j].push_back(0);
            }
            for(int k = 0; k < contain_counts[i][j].size(); k++){
                contain_counts[0][j][k] += contain_counts[i][j][k];
            }
        }
    }

	return contain_counts[0];
}
//##############################################################################
//COUNT CELL FUNCTION FILTERED

std::vector<std::vector<vertex_index_t>> filtered_count_cells(compressed_directed_graph_t& graph) {
	filtered_directed_flag_complex_t complex(graph);


   std::vector<vertex_index_t> do_vertices;
   for(int i = 0; i < graph.vertex_number(); i++){ do_vertices.push_back(i); }

    std::vector<std::vector<std::vector<vertex_index_t>>> contain_counts(PARALLEL_THREADS,
			                                                     std::vector<std::vector<vertex_index_t>>(graph.max_filtration()+1,
																														std::vector<vertex_index_t>(0)));

	std::array<cell_counter_t*, PARALLEL_THREADS> cell_counter;
	for (int i = 0; i < PARALLEL_THREADS; i++)
		cell_counter[i] = new cell_counter_t();

		complex.for_each_cell(cell_counter, do_vertices, contain_counts, 0, 10000);


  for(int i = 1; i < contain_counts.size(); i++){
        for(int j = 0; j < contain_counts[i].size(); j++){
            while(contain_counts[0][j].size() < contain_counts[i][j].size()){
                contain_counts[0][j].push_back(0);
            }
            for(int k = 0; k < contain_counts[i][j].size(); k++){
                contain_counts[0][j][k] += contain_counts[i][j][k];
            }
        }
    }

		for (int i = 0; i < PARALLEL_THREADS; i++)
		{
			std::cout<<"thread "<< i <<" size "<<cell_counter[i]->cell_count().size()<< " : " ;
			for (int j=0; j< cell_counter[i]->cell_count().size();j++)
			std::cout<< " "<< cell_counter[i]->cell_count()[j];
			std::cout<<std::endl;
		}
	return contain_counts[0];
}
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
