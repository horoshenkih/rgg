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
#include <set>
#include <unordered_set>
#include <algorithm>

using std::string;
using std::vector;
using std::pair;
using std::map;
using std::set;
using std::unordered_set;

typedef string Node;
typedef vector<double> Coordinates;

class Edge {
private:
    pair<Node, Node> node_pair;
public:
    Edge(string, string);
    bool operator==(const Edge& other) const;
    string repr() const;
};

struct hashEdge {
    std::hash<string> hasher;
    size_t operator()(const Edge& e) const {
        return hasher(e.repr());
    }
};

class NodeDescription {
private:
    unsigned int degree;
    int component_id;
    Coordinates coordinates {0., 0.};
public:
    NodeDescription();
    void increment_degree();
    unsigned int get_degree();
    int get_component_id();
    void set_component_id(int);
    Coordinates get_coordinates() const;
    void set_coordinates(const Coordinates &);
};

typedef map<Node, NodeDescription> NodeMap;
typedef map< Node, set<Node> > AdjMap;

class Nodes {
private:
    NodeMap nodes;
    typedef vector<Node> NodeContainter;
    NodeContainter nodes_list;
public:
    void add_node(const Node&);
    bool exists(const Node&) const;
    void increment_degree(Node);
    unsigned int size() const;

    NodeDescription get_description(const Node &) const;
    void set_description(const Node &, const NodeDescription &);

    typedef NodeContainter::const_iterator const_iterator;
    const_iterator begin() const {return nodes_list.begin();}
    const_iterator end() const {return nodes_list.end();}
};

class Graph {
private:
    unordered_set<Edge, hashEdge> edges;
    //NodeMap nodes;
    Nodes nodes;
    AdjMap adj_map;
public:
    const Nodes& get_nodes() const;
    unsigned int number_of_nodes() const;
    bool has_node(const Node&) const;
    unsigned int number_of_edges() const;

    void add_node(const string&);
    void add_edge(const string&, const string&);
    set<Node> neighbors(const Node&) const;
    set<Node> core_nodes(double) const;
    template <typename TNodeContainer> Graph* subgraph(const TNodeContainer &) const;
    template <typename TNodeContainer> set<Node> fringe_nodes(const TNodeContainer &) const;
    Graph* large_component() const;
    bool is_connected() const;

    NodeDescription get_node_description(const Node &) const;
    void set_node_description(const Node &, const NodeDescription &);
};

template <typename TNodeContainer>
set<Node> Graph::fringe_nodes(const TNodeContainer &nodes_container) const {
    set<Node> fringe;
    for (auto v : nodes_container) {
        for (auto n: this->neighbors(v)) {
            fringe.insert(n);
        }
    }
    for (auto node : nodes_container) {
        fringe.erase(node);
    }
    return fringe;
}

class Components {
private:
    const Graph& G;
    map<Node, int> component_ids;
    map<int, unsigned int> component_sizes;
    int components_count;
    bool is_labelled_node(const Node &) const;
    void ccR(const Node &);
public:
    Components(const Graph&);
    int node_component_id(const Node &);
    unsigned int component_size(int);
    int get_components_count();
    int max_component_id();
};