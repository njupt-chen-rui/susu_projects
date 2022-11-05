#include<iostream>
#include<set>
#include "read.h"
// #define DEBUG
// using namespace std;
// using index = std::array<int, 3>;
template<typename index>
void build_edge(vector<vector<int>> &edges, vector<vert> &nodes, vector<index> &indexes){
    vector<int> tmp;
    for(size_t i=0;i<nodes.size();i++) edges.push_back(tmp);
    for(size_t i=0;i<edges.size();i++) edges[i].push_back(-1);
    for(int i=1;i<indexes.size();i++){
        edges[indexes[i][0]].push_back(indexes[i][1]);
        edges[indexes[i][0]].push_back(indexes[i][2]);
        edges[indexes[i][1]].push_back(indexes[i][0]);
        edges[indexes[i][1]].push_back(indexes[i][2]);
        edges[indexes[i][2]].push_back(indexes[i][0]);
        edges[indexes[i][2]].push_back(indexes[i][1]);
    }
    for(size_t i=1;i<edges.size();i++){
        set<int> s(edges[i].begin(),edges[i].end());
        edges[i].assign(s.begin(),s.end());
    }
}
int main(){
    using index = std::array<int, 3>;
    std::vector<vert> P, nodes;
    std::vector<index> indexes;
    std::vector<std::vector<int>> edges;
    // read data
    string input_str;
    string ply_path = "./P.txt";
    std::cout << "deafult ply file is ./P.txt, do you want to use it? yes or no" <<endl;
    std::cin >> input_str;
    if(input_str == "yes"){
    }else if(input_str == "no"){
        std::cout << "please input the path of ply file:" << std::endl;
        std::cin >> ply_path;
    }else{
        std::cout << "please input yes or no";
        exit(1);
    }
    readply(ply_path, P);

    string inp_path = "./triangular-mesh.inp";
    std::cout << "deafult inp file is ./triangular-mesh.inp, do you want to use it? yes or no" <<endl;
    std::cin >> input_str;
    if(input_str == "yes"){
    }else if(input_str == "no"){
        std::cout << "please input the path of inp file:" << std::endl;
        std::cin >> inp_path;
    }else{
        std::cout << "please input yes or no";
        exit(1);
    } 
    readinp(inp_path, nodes, indexes);
    
    //
    build_edge(edges, nodes, indexes);
    KD_tree<index> kd_tree(P, indexes, nodes, edges);
#ifdef DEBUG
    kd_tree.build(1, (int)P.size()-1);
    kd_tree.iterative_update();
#endif
    return 0;
}