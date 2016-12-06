//
// Created by serkh on 12/6/16.
//

#ifndef HYPRG_OPTIMIZATION_H
#define HYPRG_OPTIMIZATION_H

#include "embedding_model.h"
#include "loss_function.h"
#include "pair_generator.h"

class Optimizer {
public:
    virtual void optimize_embedding(EmbeddingModel*, SmoothEdgeLoss*, PairGenerator*) = 0;
};

class SGD : public Optimizer {
private:
    double learning_rate;
    unsigned n_epoch;
    bool verbose;
public:
    SGD(double learning_rate, unsigned n_epoch, bool verbose) : learning_rate(learning_rate), n_epoch(n_epoch), verbose(verbose) {};
    void optimize_embedding(EmbeddingModel*, SmoothEdgeLoss*, PairGenerator*);
};

#endif //HYPRG_OPTIMIZATION_H
