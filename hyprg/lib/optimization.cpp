//
// Created by serkh on 12/6/16.
//

#include <iostream>
#include <cmath>
#include <map>
#include <exception>
#include <ctime>

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
    if (nc.find(node) == nc.end()) { nc[node].num_updates = 1; } else { nc[node].num_updates++; }
    //if (nc.find(node) == nc.end()) { nc[node].num_updates = 0; }
    int n = G->number_of_nodes();
    int deg = G->get_node_description(node).get_degree();
    double min_neighbors_done = std::max(1., deg/2.);
    if (nc[node].num_updates > 1000) { // TODO hardcode 5000
        nc[node].state = DONE;
    } else if (deg >= sqrt(n)) {
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
    //Node pivot;
    //int last_neighbor_degree = 0; //G->number_of_nodes();
    double max_neighbor_jaccard_similarity = 0.;
    Node max_jaccard_similar_neighbor;

    vector<Coordinates> done_neighbors;
    vector<Node> done_neighbors_nodes;
    vector<double> done_neighbors_weights;

    vector<Coordinates> processing_neighbors;

    set<Node> neighbors = G->neighbors(node);
    for (auto neighbor : neighbors) {
        int neighbor_degree = G->get_node_description(neighbor).get_degree();
        if (ns[neighbor].state == DONE) {
            done_neighbors.push_back(embedding_model->get_node_embedding(neighbor));
            done_neighbors_nodes.push_back(neighbor);
            done_neighbors_weights.push_back(1.);
            //continue;
        }
        // EXPERIMENT
            /*
        else if (ns[neighbor].state == PROCESSING and ns[neighbor].num_updates > 1000) {
            processing_neighbors.push_back(embedding_model->get_node_embedding(neighbor));
        }
             */
        /*
        if (neighbor_degree > last_neighbor_degree) {
            last_neighbor_degree = neighbor_degree;
            //pivot = neighbor;
        }
         */

        // EXPERIMENT
        /*
        set<Node> second_neighbors = G->neighbors(neighbor);
        int neighbor_common_n = 0;
        int neighbor_union_n = neighbors.size();
        for (auto sn : second_neighbors) {
            if (sn != node && neighbors.find(sn) != neighbors.end()) {
                neighbor_common_n++;
            } else {
                neighbor_union_n++;
            }
        }
        double jaccard_similarity = double(neighbor_common_n) / neighbor_union_n;//G->get_node_description(neighbor).get_degree();
        if (jaccard_similarity > max_neighbor_jaccard_similarity) {
            max_neighbor_jaccard_similarity = jaccard_similarity;
            max_jaccard_similar_neighbor = neighbor;
        }
        done_neighbors_weights.push_back(jaccard_similarity);
         */
    }
    Coordinates new_coords;
    if (!done_neighbors.size()) {
        new_coords = embedding_model->generate_random_coordinates(node);
        /*
        if (!processing_neighbors.size()) {
            new_coords = embedding_model->generate_random_coordinates(node);
        } else {
            Coordinates avg = poincare_average(processing_neighbors);
            new_coords = embedding_model->generate_random_coordinates(node);
            new_coords[1] = avg[1];
        }
         */
    } else {
        //new_coords = embedding_model->generate_random_coordinates(node, pivot);
        /*
        for (auto dn : done_neighbors_nodes) {
            int deg = G->get_node_description(dn).get_degree();
            done_neighbors_weights.push_back(1./ deg);
        }
         */
        Coordinates avg = poincare_average(done_neighbors, done_neighbors_weights);
        //Coordinates avg = poincare_average(done_neighbors);
        new_coords = embedding_model->generate_random_coordinates(node);
        new_coords[1] = avg[1];
    }
    if (new_coords[0] < 0) {
        throw std::domain_error("(in optimization, __generate_new_cordinates_pivot) generated point has negative radius!");
    }
    if (std::isnan(new_coords[1])) {
        throw std::domain_error("(in optimization, __generate_new_cordinates_pivot) generated point angle is NaN!");
    }
    return new_coords;
}

void SGD::optimize_embedding(
        PoincareModel *embedding_model,
        SmoothEdgeLoss *loss_function,
        PairGenerator *pair_generator,
        Graph* G
) {
    int total_iterations=0;
    NodeStates node_update_count;
    unsigned n = embedding_model->get_node_count();
    double n_updates_shift = 10;
    double gradient_step_time_exponent = 0.75;
    double gradient_step_degree_exponent = 0.5;

    // EXPERIMENT
    double max_beta = 5; //1 * log(n);
    double min_beta = 5; //1 * log(n);
    loss_function->set_beta(max_beta);
    if (verbose) {
        std::cout << "n_epoch=" << n_epoch << std::endl;
        std::cout << "Generate pairs" << std::endl;
    }
    pair_generator->generate_pairs();
    for (int epoch=0; epoch < n_epoch; ++epoch) {
        clock_t epoch_begin = std::clock();
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
            if (!embedding_model->is_external_init()) {
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

                // TODO EXPERIMENT
                if (node_update_count[nodepair.first].state == PENDING or node_update_count[nodepair.second].state == PENDING) {
                    continue;
                }
                // /EXPERIMENT
            }
            DistanceInfo di = embedding_model->get_edge_distance_info(edge);
            double loss = loss_function->loss(edge, di);
            total_loss += loss;
            //if (!embedding_model->is_node_embedded(nodepair.first) or !embedding_model->is_node_embedded(nodepair.second)) {
            //    continue;
            //}

            auto grad = loss_function->loss_gradient(edge, di);
            auto x = embedding_model->get_edge_vector(edge);

            double first_multiplier = 1. / pow(node_update_count[nodepair.first].num_updates+n_updates_shift, gradient_step_time_exponent);
            double second_multiplier = 1. / pow(node_update_count[nodepair.second].num_updates+n_updates_shift, gradient_step_time_exponent);
            int first_degree = G->get_node_description(nodepair.first).get_degree();
            int second_degree = G->get_node_description(nodepair.second).get_degree();
            first_multiplier /= pow(first_degree, gradient_step_degree_exponent);
            second_multiplier /= pow(second_degree, gradient_step_degree_exponent);
            // Don't touch core
            if (first_degree >= sqrt(n)) {
                first_multiplier = 0.;
            }
            if (second_degree >= sqrt(n)) {
                second_multiplier = 0.;
            }
            for (int i = 0; i < x.size(); ++i) {
                x[i] -= learning_rate * grad[i] * (i < 2 ? first_multiplier : second_multiplier) / pow(++total_iterations, 0.25);
            }

            // TODO hack for nodes with degree=1
            /*
            if (G->get_node_description(nodepair.first).get_degree() == 1) {
                auto neighs = G->neighbors(nodepair.first);
                for (auto neigh : neighs) {
                    Coordinates neigh_embedding = embedding_model->get_node_embedding(neigh);
                    x[1] = neigh_embedding[1];
                }
            }
            if (G->get_node_description(nodepair.second).get_degree() == 1) {
                auto neighs = G->neighbors(nodepair.second);
                for (auto neigh : neighs) {
                    Coordinates neigh_embedding = embedding_model->get_node_embedding(neigh);
                    x[3] = neigh_embedding[1];
                }
            }
             */
            // /TODO hack
            //if (embedding_model->is_external_init() or node_update_count[nodepair.first].state == PROCESSING) {
            if (node_update_count[nodepair.first].state == PROCESSING) {
                embedding_model->set_node_embedding(nodepair.first,  Coordinates{std::max(0., x[0]), x[1]});
                node_update_count[nodepair.first].num_updates++;
            }
            //if (embedding_model->is_external_init() or node_update_count[nodepair.second].state == PROCESSING) {
            if (node_update_count[nodepair.second].state == PROCESSING) {
                embedding_model->set_node_embedding(nodepair.second, Coordinates{std::max(0., x[2]), x[3]});
                node_update_count[nodepair.second].num_updates++;
            }
        }
        if (verbose) {
            std::cout << "Average loss: " << total_loss / N << std::endl;
            std::cout << "Elapsed time: " << float(std::clock() - epoch_begin) / CLOCKS_PER_SEC << std::endl;
            std::cout << "===========" << std::endl;
            if (true) {
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
        }
        double beta = loss_function->get_beta();
        loss_function->set_beta(min_beta + 0.9 *(beta-min_beta));
    }
}