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
//#include <boost/graph/adjacency_list.hpp>

using std::string;
using std::vector;
using std::pair;
using std::map;
using std::set;

typedef string Node;

class Edge {
private:
    pair<Node, Node> node_pair;
public:
    Edge(string, string);
};

class NodeDescription {
private:
    unsigned int degree;
    int component_id;
public:
    NodeDescription();
    void increment_degree();
    int get_component_id();
    void set_component_id(int);
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
    bool exists(const Node&);
    void increment_degree(Node);
    unsigned int size() const;

    typedef NodeContainter::const_iterator const_iterator;
    const_iterator begin() const {return nodes_list.begin();}
    const_iterator end() const {return nodes_list.end();}
};

class Graph {
private:
    vector<Edge> edges;
    //NodeMap nodes;
    Nodes nodes;
    AdjMap adj_map;
public:
    const Nodes& get_nodes() const;
    unsigned int number_of_nodes() const;
    unsigned int number_of_edges() const;
    void add_node(const string&);
    void add_edge(const string&, const string&);
    set<Node> neighbors(const Node&) const;
    //set<Node> neighbors(const string&);
};

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
};