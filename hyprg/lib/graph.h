//
// Created by serkh on 11/12/16.
//

#ifndef HYPRG_GRAPH_H
#define HYPRG_GRAPH_H

#endif //HYPRG_GRAPH_H
#include <string>
#include <vector>
#include <utility>
#include <map>

using std::string;
using std::vector;
using std::pair;
using std::map;

typedef string Node;
typedef map<Node, unsigned int> CounterMap;

class Edge {
private:
    pair<Node, Node> node_pair;
public:
    Edge(string, string);
};

class Graph {
private:
    vector<Edge> edges;
    CounterMap nodes;
public:
    unsigned int number_of_nodes();
    unsigned int number_of_edges();
    void add_node(string);
    void add_edge(string, string);
};