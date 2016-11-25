//
// Created by serkh on 11/12/16.
//

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "graph.h"

using std::string;

Graph read_graph_from_file(const char *filename) {
    Graph G;
    std::ifstream infile(filename);
    string line;
    while (std::getline(infile, line)) {
        if (line.find("#") == 0) {
            continue;
        }
        std::istringstream iss(line);
        string a, b;
        if (!(iss >> a >> b)) { throw std::runtime_error("wrong line format"); }
        G.add_edge(a, b);
    }
    return G;
}

