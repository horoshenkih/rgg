//
// Created by serkh on 12/6/16.
//

#ifndef HYPRG_EMBEDDING_MODEL_H
#define HYPRG_EMBEDDING_MODEL_H

#include <map>
#include <vector>
#include <random>

#include "graph.h"
#include "pair_generator.h"

typedef std::vector<double> Gradient;
typedef std::map<Node, Gradient> DistanceGradients;

class DistanceInfo {
private:
    double radius;
    double distance;
    DistanceGradients distance_gradients;
public:
    DistanceInfo(double r, double d, DistanceGradients dg) : radius(r), distance(d), distance_gradients(dg) {}
    double get_radius() const { return radius; }
    double get_distance() const { return distance; }
    DistanceGradients get_distance_gradients() const { return distance_gradients; }
};

// TODO slow...
/*
class RandomUniformGenerator {
private:
    double a, b;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
public:
    RandomUniformGenerator(double a, double b) : a(a), b(b) {
        std::random_device rd;
        gen = std::mt19937(rd());
        dis = std::uniform_real_distribution<>(a,b);
    }
    double operator()() { return dis(gen); }
};
*/

typedef std::map<Node, Coordinates> NodeEmbeddings;
class EmbeddingModel {
protected:
    NodeEmbeddings node_coordinates;
public:
    bool is_node_embedded(const Node &) const;
    Coordinates get_node_embedding(const Node &) const;
    void set_node_embedding(const Node &, Coordinates);
    Coordinates get_edge_vector(const Edge &) const;
    NodeEmbeddings get_all_node_embeddings() const;
    unsigned get_node_count() const;
    virtual Coordinates get_state_vector() const = 0;  // TODO need it?
    virtual void set_state_vector(const Coordinates &) = 0; // TODO need it?
    virtual DistanceInfo get_edge_distance_info(const Edge &) const = 0;
    virtual Coordinates generate_random_coordinates(const Node& node) const = 0;
    virtual Coordinates generate_random_coordinates(const Node& node, const Node& pivot) const = 0;
    virtual Coordinates generate_random_coordinates(const Node& node, const Coordinates& pivot_coordinates) const = 0;
};

class PoincareModel : public EmbeddingModel {
private:
    Graph* G;
    double radius;
public:
    PoincareModel(Graph* G);
    void set_state_vector(const Coordinates &);
    Coordinates get_state_vector() const;
    DistanceInfo get_edge_distance_info(const Edge &) const;
    double generate_radius(const Node&) const;
    Coordinates generate_random_coordinates(const Node& node) const;
    Coordinates generate_random_coordinates(const Node&, const Node& pivot) const;
    Coordinates generate_random_coordinates(const Node&, const Coordinates& pivot_coordinates) const;
};
#endif //HYPRG_EMBEDDING_MODEL_H
