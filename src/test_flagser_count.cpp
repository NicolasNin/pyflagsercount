#include "filtered_flagser-count.cpp"

int main(){

// std::vector<vertex_index_t> row   {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4};
// std::vector<vertex_index_t> col{1, 3, 4, 2, 3, 4, 0, 1, 3, 4, 2, 4, 0, 2, 3};
// int N=5;
std::vector<vertex_index_t> row{0,0,1};
std::vector<vertex_index_t> col{1,2,2};
int N=3;
auto graph = directed_graph_t(N);
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
auto a=grow_trees(graph);
int c=0;
while(c<10){
    c++;
for (int i=0;i<10;i++)
    a=    grow_trees(graph);
}
    return 0;
}