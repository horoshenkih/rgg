//
// Created by serkh on 12/5/16.
//

#ifndef HYPRG_PAIR_GENERATOR_H
#define HYPRG_PAIR_GENERATOR_H

#include <vector>
#include "graph.h"

typedef vector<Edge> WeightedPairs;
class PairGenerator {
private:
    Graph* G;
    WeightedPairs pairs;
public:
    PairGenerator(Graph* G);
    const WeightedPairs& get_pairs() const;
    void generate_pairs(
            double ratio_to_second=2.,
            double ratio_between_first=1.,
            double ratio_first_second=1.,
            double ratio_between_second=1.,
            double ratio_random=1.
    );
    void shuffle_pairs();
};


#endif //HYPRG_PAIR_GENERATOR_H
