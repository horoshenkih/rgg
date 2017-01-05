//
// Created by serkh on 12/6/16.
//

#include <cmath>
#include <algorithm>

#include "embedding_model.h"
#include "hyperbolic.h"
#include "utils.h"

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

unsigned EmbeddingModel::get_node_count() const { return node_coordinates.size(); }

Coordinates EmbeddingModel::get_edge_vector(const Edge &edge) const {
    pair<Node, Node> node_pair = edge.get_node_pair();
    Coordinates v1 = this->get_node_embedding(node_pair.first);
    Coordinates v2 = this->get_node_embedding(node_pair.second);
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

PoincareModel::PoincareModel(Graph* G) {
    this->G = G;
    unsigned n = G->number_of_nodes();
    radius = 2 * log(n);

    for (auto v : G->get_nodes()) {
        /*
        NodeDescription d = G->get_node_description(v);
        unsigned int deg = d.get_degree();
        double r = 2. * log(double(n) / deg);
        double r_jitter = 2. / sqrt(deg) * (2. * static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 1.);
        r += r_jitter;
        float phi = 2 * 3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        Coordinates c{r, phi};
        node_coordinates[v] = c;
         */
        node_coordinates[v] = generate_random_coordinates(v);
    }
    external_init = false;
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


double PoincareModel::generate_radius(const Node &node) const {
    unsigned n = G->number_of_nodes();
    NodeDescription d = G->get_node_description(node);
    unsigned int deg = d.get_degree();
    //double r = generate_random_uniform(0, 2. * log(double(n) / deg));
    double r = 2. * log(double(n) / deg);
    double r_jitter = std::min(2./sqrt(deg), r - 1./n);
    r += generate_random_uniform(-r_jitter, r_jitter);
    return std::max(0., std::min(r, radius));
}

Coordinates PoincareModel::generate_random_coordinates(const Node &node) const {
    double r = generate_radius(node);
    double phi = generate_random_uniform(0, 2*M_PI);
    return Coordinates{r, phi};
}

Coordinates PoincareModel::generate_random_coordinates(const Node &node, const Coordinates &pivot_coordinates) const {
    //Coordinates pivot_coordinates = get_node_embedding(pivot);
    double r_n = pivot_coordinates[0];
    double phi_n = pivot_coordinates[1];
    double cos_delta_phi = (pow(cosh(r_n), 2) - cosh(radius)) / pow(sinh(r_n), 2);
    cos_delta_phi = std::max(-0.999, std::min(0.999, cos_delta_phi));
    double delta_phi = acos(cos_delta_phi);

    double r = generate_radius(node);
    double phi = generate_random_uniform(phi_n - delta_phi, phi_n + delta_phi);
    return Coordinates{r, phi};
}

Coordinates PoincareModel::generate_random_coordinates(const Node &node, const Node &pivot) const {
    Coordinates pivot_coordinates = get_node_embedding(pivot);
    return generate_random_coordinates(node, pivot_coordinates);
    /*
    double r_n = pivot_coordinates[0];
    double phi_n = pivot_coordinates[1];
    double cos_delta_phi = (pow(cosh(r_n), 2) - cosh(radius)) / pow(sinh(r_n), 2);
    cos_delta_phi = std::max(-0.999, std::min(0.999, cos_delta_phi));
    double delta_phi = acos(cos_delta_phi);

    double r = generate_radius(node);
    double phi = generate_random_uniform(phi_n - delta_phi, phi_n + delta_phi);
    return Coordinates{r, phi};
     */
}

void PoincareModel::set_external_init(bool ei) { external_init = ei; }
bool PoincareModel::is_external_init() {return external_init; }

Coordinates PoincareModel::get_state_vector() const {
    // TODO
    return Coordinates{0,0};
}

void PoincareModel::set_state_vector(const Coordinates &) {
    // TODO
}
