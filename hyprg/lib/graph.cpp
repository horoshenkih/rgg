//
// Created by serkh on 11/12/16.
//
#include "graph.h"

Edge::Edge(string a, string b) {
    node_pair = pair<Node, Node>(Node(a), Node(b));
}

unsigned int Graph::number_of_nodes() {
    return nodes.size();
}
unsigned int Graph::number_of_edges() {
    return edges.size();
}

void Graph::add_node(string n) {
    Node node(n);// = Node(n);
    if (nodes.find(node) != nodes.end()) {
        nodes[node] += 1;
    } else {
        nodes[node] = 1;
    }
}

void Graph::add_edge(string a, string b) {
    edges.push_back(Edge(a, b));
    this->add_node(a);
    this->add_node(b);
}
