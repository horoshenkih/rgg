//
// Created by serkh on 12/6/16.
//

#include <iostream>
#include <cmath>
#include <map>
#include <exception>
#include "optimization.h"
#include "hyperbolic.h"

enum State { PENDING, PROCESSING, DONE};

// TODO combine NodeStates and "private" functions to class
struct NodeOptimizationState {
    State state;
    int num_updates;
};

typedef std::map<Node, NodeOptimizationState> NodeStates;
void __update_node_state(NodeStates &nc, const Node & node, Graph* G) {
    //if (nc.find(node) == nc.end()) { nc[node].num_updates = 1; } else { nc[node].num_updates++; }
    if (nc.find(node) == nc.end()) { nc[node].num_updates = 0; }
    int n = G->number_of_nodes();
    int deg = G->get_node_description(node).get_degree();
    double min_neighbors_done = std::max(1., deg/2.);
    if (nc[node].num_updates > 5000) { // TODO hardcode
        nc[node].state = DONE;
    } else if (deg > sqrt(n)) {
        nc[node].state = PROCESSING;
    } else {
        bool is_any_neighbor_done = false;
        //int n_neighbors_done = 0;
        for (auto neigh : G->neighbors(node)) {
            if (nc[neigh].state == DONE) {
                is_any_neighbor_done = true;
                //n_neighbors_done++;
                break;
            }
        }
        if (is_any_neighbor_done) {
        //if (n_neighbors_done >= min_neighbors_done) {
            nc[node].state = PROCESSING;
        } else {
            nc[node].state = PENDING;
        }
    }
}

Coordinates __generate_new_coordinates_pivot(const Node& node, NodeStates& ns, EmbeddingModel* embedding_model, Graph* G) {
    Node pivot;
    int last_neighbor_degree = 0; //G->number_of_nodes();
    vector<Coordinates> done_neighbors;
    vector<Node> done_neighbors_nodes;
    for (auto neighbor : G->neighbors(node)) {
        if (ns[neighbor].state != DONE) {
            continue;
        }
        int neighbor_degree = G->get_node_description(neighbor).get_degree();
        if (neighbor_degree > last_neighbor_degree) {
            last_neighbor_degree = neighbor_degree;
            pivot = neighbor;
        }
        done_neighbors.push_back(embedding_model->get_node_embedding(neighbor));
        done_neighbors_nodes.push_back(neighbor);
    }
    Coordinates new_coords;
    if (last_neighbor_degree == 0) {
        new_coords = embedding_model->generate_random_coordinates(node);
    } else {
        //new_coords = embedding_model->generate_random_coordinates(node, pivot);
        /*
        vector<double> done_neighbors_weights;
        for (auto dn : done_neighbors_nodes) {
            int deg = G->get_node_description(dn).get_degree();
            done_neighbors_weights.push_back(1./ deg);
        }
        Coordinates avg = poincare_average(done_neighbors, done_neighbors_weights);
         */
        Coordinates avg = poincare_average(done_neighbors);
        new_coords = embedding_model->generate_random_coordinates(node);
        new_coords[1] = avg[1];
    }
    return new_coords;
}

void SGD::optimize_embedding(EmbeddingModel *embedding_model, SmoothEdgeLoss *loss_function, PairGenerator *pair_generator, Graph* G) {
    int total_iterations=0;
    NodeStates node_update_count;
    unsigned n = embedding_model->get_node_count();
    double gradient_step_time_exponent = 0.75;

    // EXPERIMENT
    double max_beta = 4 * log(n);
    double min_beta = 1 * log(n);
    loss_function->set_beta(max_beta);
    pair_generator->generate_pairs(30,30,0);
    for (int epoch=0; epoch < n_epoch; ++epoch) {
        pair_generator->shuffle_pairs();
        if (verbose) {
            std::cout << "Epoch " << epoch+1 << " / " << n_epoch << std::endl;
            std::cout << "beta=" << loss_function->get_beta() << std::endl;
        }
        double total_loss = 0;
        unsigned N = 0;
        for (auto edge: pair_generator->get_pairs()) {
            N++;
            pair<Node, Node> nodepair = edge.get_node_pair();
            // EXPERIMENT
            __update_node_state(node_update_count, nodepair.first, G);
            __update_node_state(node_update_count, nodepair.second, G);
            if (node_update_count[nodepair.first].state == PROCESSING and node_update_count[nodepair.first].num_updates == 0) {
                Coordinates new_c = __generate_new_coordinates_pivot(nodepair.first, node_update_count, embedding_model, G);
                embedding_model->set_node_embedding(nodepair.first, new_c);
            }
            if (node_update_count[nodepair.second].state == PROCESSING and node_update_count[nodepair.second].num_updates == 0) {
                Coordinates new_c = __generate_new_coordinates_pivot(nodepair.second, node_update_count, embedding_model, G);
                embedding_model->set_node_embedding(nodepair.second, new_c);
            }
            //if (node_update_count[nodepair.first].state != PROCESSING and node_update_count[nodepair.second].state != PROCESSING) {
            //    continue;
            //}
            // /EXPERIMENT
            DistanceInfo di = embedding_model->get_edge_distance_info(edge);
            double loss = loss_function->loss(edge, di);
            total_loss += loss;
            //if (!embedding_model->is_node_embedded(nodepair.first) or !embedding_model->is_node_embedded(nodepair.second)) {
            //    continue;
            //}

            auto grad = loss_function->loss_gradient(edge, di);
            auto x = embedding_model->get_edge_vector(edge);

            double first_multiplier = 1. / pow(node_update_count[nodepair.first].num_updates+1, gradient_step_time_exponent);
            double second_multiplier = 1. / pow(node_update_count[nodepair.second].num_updates+1, gradient_step_time_exponent);
            for (int i = 0; i < x.size(); ++i) {
                x[i] -= learning_rate * grad[i] * (i < 2 ? first_multiplier : second_multiplier); // / pow(++total_iterations, 0.25);
            }

            // EXPERIMENT
            if (node_update_count[nodepair.first].state == PROCESSING) {
                embedding_model->set_node_embedding(nodepair.first,  Coordinates{std::max(0., x[0]), x[1]});
                node_update_count[nodepair.first].num_updates++;
            }
            if (node_update_count[nodepair.second].state == PROCESSING) {
                embedding_model->set_node_embedding(nodepair.second, Coordinates{std::max(0., x[2]), x[3]});
                node_update_count[nodepair.second].num_updates++;
            }
        }
        if (verbose) {
            std::cout << "Average loss: " << total_loss / N << std::endl;
            unsigned n_pending = 0;
            unsigned n_processing = 0;
            unsigned n_done = 0;
            for (auto nc : node_update_count) {
                switch (nc.second.state) {
                    case PENDING: n_pending++; break;
                    case PROCESSING: n_processing++; break;
                    case DONE: n_done++; break;
                }
            }
            std::cout << "Total nodes: " << node_update_count.size() << std::endl;
            std::cout << "pending: " << n_pending << std::endl;
            std::cout << "processing: " << n_processing << std::endl;
            std::cout << "done: " << n_done << std::endl;
        }
        double beta = loss_function->get_beta();
        loss_function->set_beta(min_beta + 0.9 *(beta-min_beta));
    }
}