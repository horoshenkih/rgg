//
// Created by serkh on 12/6/16.
//

#include <cmath>

#include "embedding_model.h"
#include "hyperbolic.h"

NodeEmbeddings EmbeddingModel::get_all_node_embeddings() const {
    return node_coordinates;
}

bool EmbeddingModel::is_node_embedded(const Node &n) const {
    return node_coordinates.find(n) != node_coordinates.end();
}

Coordinates EmbeddingModel::get_node_embedding(const Node &n) const {
    if (is_node_embedded(n)) {
        return node_coordinates.find(n)->second;
    } else {
        return Coordinates{0,0};
    }
}

void EmbeddingModel::set_node_embedding(const Node &node, Coordinates coords) {
    node_coordinates[node] = coords;
}

Coordinates EmbeddingModel::get_edge_vector(const Edge &edge) const {
    pair<Node, Node> node_pair = edge.get_node_pair();
    Coordinates v1 = this->get_node_embedding(node_pair.first);
    Coordinates v2 = this->get_node_embedding(node_pair.second);
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

PoincareModel::PoincareModel(Graph* G) {
    unsigned n = G->number_of_nodes();
    radius = 2 * log(n);

    for (auto v : G->get_nodes()) {
        NodeDescription d = G->get_node_description(v);
        unsigned int deg = d.get_degree();
        double r = 2. * log(double(n) / deg);
        float phi = 2 * 3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        Coordinates c{r, phi};
        node_coordinates[v] = c;
    }
}

DistanceInfo PoincareModel::get_edge_distance_info(const Edge &edge) const {
    auto node_pair = edge.get_node_pair();
    Coordinates c1 = this->get_node_embedding(node_pair.first);
    Coordinates c2 = this->get_node_embedding(node_pair.second);

    double d = distance(c1[0], c1[1], c2[0], c2[1]);
    vector<double> grad_d = grad_distance(c1[0], c1[1], c2[0], c2[1]);

    DistanceGradients dg;
    dg[node_pair.first]  = vector<double>{grad_d[0], grad_d[1]};
    dg[node_pair.second] = vector<double>{grad_d[2], grad_d[3]};
    return DistanceInfo(radius, d, dg);
}

Coordinates PoincareModel::get_state_vector() const {
    // TODO
    return Coordinates{0,0};
}

void PoincareModel::set_state_vector(const Coordinates &) {
    // TODO
}