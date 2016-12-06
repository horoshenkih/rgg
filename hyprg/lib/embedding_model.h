//
// Created by serkh on 12/6/16.
//

#ifndef HYPRG_EMBEDDING_MODEL_H
#define HYPRG_EMBEDDING_MODEL_H

#include <map>
#include "graph.h"
#include "pair_generator.h"

class DistanceInfo {

};

typedef std::map<Node, Coordinates> NodeEmbeddings;
class EmbeddingModel {
protected:
     NodeEmbeddings node_coordinates;
public:
    bool is_node_embedded(const Node &) const;
    Coordinates get_node_embedding(const Node &) const;
    NodeEmbeddings get_all_node_embeddings() const;
    virtual Coordinates get_state_vector() const = 0;
    virtual void set_state_vector(const Coordinates &) = 0;
    //virtual DistanceInfo get_distance_info(const WeightedPairs &) = 0;
};

class PoincareModel : public EmbeddingModel {
public:
    PoincareModel(Graph &G);
    void set_state_vector(const Coordinates &);
    Coordinates get_state_vector() const;
};
#endif //HYPRG_EMBEDDING_MODEL_H
