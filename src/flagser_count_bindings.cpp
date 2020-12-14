#include <stdio.h>
#include <iostream>

#include "filtered_flagser-count.cpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(pyflagsercontain2, m) {

  m.doc() = "Python interface for flagser_count";

  m.def("compute_cell_count", [](vertex_index_t num_vertices,
                                 std::vector<std::vector<value_t>>& edges) {
    // Save std::cout status
    auto cout_buff = std::cout.rdbuf();

    // Building the filtered directed graph
    auto graph = directed_graph_t(num_vertices);

    for (auto& edge : edges) {
        graph.add_edge(edge[0], edge[1]);
    }

    // Disable cout
    std::cout.rdbuf(nullptr);

    // Running flagser-count's count_cells routine
    auto cell_count = count_cells(graph);

    // Re-enable again cout
    std::cout.rdbuf(cout_buff);

    return cell_count;
  });

    m.def("compute_cell_count_filtered", [](vertex_index_t num_vertices,
                                 std::vector<std::vector<value_t>>& edges) {
    // Save std::cout status
    auto cout_buff = std::cout.rdbuf();

    // Building the filtered directed graph
    auto graph = compressed_directed_graph_t(num_vertices);

    for (auto& edge : edges) {
        graph.add_edge_weighted(edge[0], edge[1],edge[2]);
    }

    // Disable cout
    std::cout.rdbuf(nullptr);

    // Running flagser-count's count_cells routine
    auto cell_count = filtered_count_cells(graph);

    // Re-enable again cout
    std::cout.rdbuf(cout_buff);

    return cell_count;
  });

}
