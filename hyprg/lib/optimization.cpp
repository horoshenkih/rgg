//
// Created by serkh on 12/6/16.
//

#include <iostream>
#include "optimization.h"

void SGD::optimize_embedding(EmbeddingModel *embedding_model, SmoothEdgeLoss *loss_function, PairGenerator *pair_generator) {
    for (int epoch=0; epoch < n_epoch; ++epoch) {
        pair_generator->shuffle_pairs();
        if (verbose) {
            std::cout << "Epoch " << epoch+1 << " / " << n_epoch << std::endl;
        }
        for (auto edge: pair_generator->get_pairs()) {
            pair<Node, Node> nodepair = edge.get_node_pair();
            if (!embedding_model->is_node_embedded(nodepair.first) or !embedding_model->is_node_embedded(nodepair.second)) {
                continue;
            }
            DistanceInfo di = embedding_model->get_edge_distance_info(edge);
            auto grad = loss_function->loss_gradient(edge, di);
            auto x = embedding_model->get_edge_vector(edge);

            for (int i = 0; i < x.size(); ++i) {
                x[i] -= learning_rate * grad[i];
            }

            embedding_model->set_node_embedding(nodepair.first,  Coordinates{x[0], x[1]});
            embedding_model->set_node_embedding(nodepair.second, Coordinates{x[2], x[3]});
        }
    }
}