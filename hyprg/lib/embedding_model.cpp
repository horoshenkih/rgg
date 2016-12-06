//
// Created by serkh on 12/6/16.
//

#include <cmath>

#include "embedding_model.h"

NodeEmbeddings EmbeddingModel::get_all_node_embeddings() const {
    return node_coordinates;
}

bool EmbeddingModel::is_node_embedded(const Node &n) const {
    return node_coordinates.find(n) != node_coordinates.end();
}

Coordinates EmbeddingModel::get_node_embedding(const Node &n) const {
    return node_coordinates.find(n)->second;
}

PoincareModel::PoincareModel(Graph& G) {
    unsigned n = G.number_of_nodes();
    for (auto v : G.get_nodes()) {
        NodeDescription d = G.get_node_description(v);
        unsigned int deg = d.get_degree();
        double r = 2. * log(double(n) / deg);
        float phi = 2 * 3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        Coordinates c{r, phi};
        node_coordinates[v] = c;
    }
}

Coordinates PoincareModel::get_state_vector() const {
    // TODO
    return Coordinates{0,0};
}

void PoincareModel::set_state_vector(const Coordinates &) {
    // TODO
}