//
// Created by serkh on 11/12/16.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "graph.h"
#include "embedding_model.h"

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

void write_embedding_to_file(const EmbeddingModel& embedding, const char *filename) {
    std::ofstream out_embeddings_file;
    out_embeddings_file.open(filename);
    for (auto node_embedding: embedding.get_all_node_embeddings()) {
        Node n = node_embedding.first;
        Coordinates coords = node_embedding.second;
        out_embeddings_file << n;
        for (double x : coords) {
            out_embeddings_file << "\t" << x;
        }
        out_embeddings_file << std::endl;
    }
    out_embeddings_file.close();
}
