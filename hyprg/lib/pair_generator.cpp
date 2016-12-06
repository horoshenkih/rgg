//
// Created by serkh on 12/5/16.
//
#include <iostream>
#include "pair_generator.h"
#include "utils.h"

PairGenerator::PairGenerator(Graph* G, double ratio_to_second, double ratio_between_first, double ratio_random) {
    // edges
    Edges current_pairs;
    for (Edge e : G->get_edges()) {
        e.set_value(1.);
        e.set_weight(1.);
        //pairs.push_back(e);
        current_pairs.insert(e);
    }
    // nedges
    int total_nedges = 0;
    auto sorted_nodes = G->get_sorted_nodes();
    for (Node n : sorted_nodes) {
        int degree = G->get_node_description(n).get_degree();
        // 1. to second
        auto first_neigh = G->neighbors(n);
        vector<Node> second_neigh;
        for (Node neigh1 : first_neigh) {
            for (Node neigh2 : G->neighbors(neigh1)) {
                if (neigh2 != n) {
                    second_neigh.push_back(neigh2);
                }
            }
        }
        shuffle_vector(second_neigh);
        int n_second_vertex_nedges = 0;
        for (Node second_n : second_neigh) {
            if (n_second_vertex_nedges > degree * ratio_between_first) {
                break;
            }
            Edge second_e(n, second_n);
            second_e.set_value(0.);

            if (current_pairs.find(second_e) == current_pairs.end()) {
                current_pairs.insert(second_e);
                n_second_vertex_nedges++;
                total_nedges++;
            }
        }

        // 2. between first
        int n_between_first_nedges = 0;
        vector<Node> first_neigh_vec(first_neigh.begin(), first_neigh.end());
        unsigned n_first_neigh = first_neigh_vec.size();
        for (int i = 0; i < n_first_neigh; ++i) {
            for (int j = 0; j < i; ++j) {
                if (n_between_first_nedges > degree * ratio_to_second) {
                    break;
                }
                Node v1 = first_neigh_vec[i];
                Node v2 = first_neigh_vec[j];
                Edge first_e(v1, v2);
                first_e.set_value(0.);

                if (current_pairs.find(first_e) == current_pairs.end()) {
                    current_pairs.insert(first_e);
                    n_between_first_nedges++;
                    total_nedges++;
                }
            }
        }

        // 3.random
        int n_random_vertices = 0;
        int max_n_random_vertices = static_cast<int>(ratio_random * degree);
        vector<Node> random_vertices;
        random_vertices.resize(max_n_random_vertices);
        auto all_nodes = G->get_nodes();
        n_random_elements(all_nodes.begin(), all_nodes.end(), random_vertices.begin(), max_n_random_vertices);
        for (auto rand_n : random_vertices) {
            if (rand_n == n) {
                continue;
            }
            if (n_random_vertices > max_n_random_vertices) {
                break;
            }

            Edge random_e(rand_n, n);
            if (current_pairs.find(random_e) == current_pairs.end()) {
                current_pairs.insert(random_e);
                n_random_vertices++;
                total_nedges++;
            }
        }
    }
    double non_edge_weight = 1.;
    if (total_nedges > 0) {
        non_edge_weight = static_cast<double>(G->number_of_edges()) / total_nedges;
    }
    // copy
    pairs = vector<Edge>(current_pairs.begin(), current_pairs.end());
    // set weights
    for (auto e_it = pairs.begin(); e_it != pairs.end(); ++e_it) {
        if (e_it->get_value() < 1) {
            e_it->set_weight(non_edge_weight);
        }
    }
    shuffle_vector(pairs);
}

const WeightedPairs& PairGenerator::get_pairs() const { return pairs; }
void PairGenerator::shuffle_pairs() { shuffle_vector(pairs); }