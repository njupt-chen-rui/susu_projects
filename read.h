#include<fstream>
#include<iostream>
#include "datastruct.h"

using namespace std;
// using index = std::array<int, 3>;
void readply(string Path, vector<vert> &P){
    ifstream fin;
    fin.open(Path);
    if (!fin.good()) {
        cout << "file open failed";
        exit(1);
    }
    double a, b, c;
    vert tmp{0, 0, 0};
    P.push_back(tmp);
    while(fin>>a){
        fin>>b>>c;
        tmp.x = a;
        tmp.y = b;
        tmp.z = c;
        P.push_back(tmp);
    }
    fin.close();
}

template<typename index>
void readinp(string Path, vector<vert> &nodes, vector<index> &indexes){
    ifstream fin;
    fin.open(Path);
    if (!fin.good()) {
        cout << "file open failed";
        exit(1);
    }
    string str, c;
    char ch;
    int a;
    double b1, b2, b3;
    int cnt = 1;
    vert tmp{0, 0, 0};
    nodes.push_back(tmp);
    while(fin>>str){
        if(str=="*NODE") break;
    }
    int ant = 0;
    while(1){
        fin.get();
        fin.get();
        ch = fin.peek();
        if(ch == '*'){
            break;
        }
        fin>>a>>c>>b1>>c>>b2>>c>>b3;
        tmp.x = b1;
        tmp.y = b2;
        tmp.z = b3;
        nodes.push_back(tmp);
    }
    while(fin>>str){
        if(str=="*ELEMENT,TYPE=S3,ELSET=mesh") break;
    }
    int meshid;
    index tmp2{0, 0, 0};
    indexes.push_back(tmp2);
    while(1){
        fin.get();
        fin.get();
        ch = fin.peek();
        if(ch == '*') {
            break;
        }
        fin>>meshid>>c;
        fin>>tmp2[0]>>c;
        fin>>tmp2[1]>>c;
        fin>>tmp2[2];
        indexes.push_back(tmp2);
    }
    fin.close();
}