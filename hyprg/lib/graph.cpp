//
// Created by serkh on 11/12/16.
//
#include "graph.h"
#include <iostream>
#include <stack>

using std::cout;
using std::endl;
using std::stack;

Edge::Edge(string a, string b) {
    Node na(a);
    Node nb(b);
    if (na < nb) {
        node_pair = pair<Node, Node>(na, nb);
    } else {
        node_pair = pair<Node, Node>(nb, na);
    }
}

bool Edge::operator==(const Edge &other) const {
    return this->node_pair.first == other.node_pair.first && this->node_pair.second == other.node_pair.second;
}

string Edge::repr() const {
    return this->node_pair.first + "\t" + this->node_pair.second;
}

bool Nodes::exists(const Node& node) const {
    return nodes.find(node) != nodes.end();
}
void Nodes::add_node(const Node& node) {
    if (!exists(node)) {
        nodes[node] = NodeDescription();
        nodes_list.push_back(node);
    }
}
unsigned int Nodes::size() const {return nodes.size();}
void Nodes::increment_degree(Node n) {
    if (exists(n)){
        nodes[n].increment_degree();
    }
}

NodeDescription::NodeDescription() : degree(0), component_id(-1) {}
void NodeDescription::increment_degree() {degree++;}
int NodeDescription::get_component_id() {return component_id;}
void NodeDescription::set_component_id(int id) {component_id = id;}

const Nodes& Graph::get_nodes() const {return nodes;}

unsigned int Graph::number_of_nodes() const {
    return nodes.size();
}
unsigned int Graph::number_of_edges() const {
    return edges.size();
}

void Graph::add_node(const string& n) {
    Node node(n);
    nodes.add_node(node);
}

void Graph::add_edge(const string& a, const string& b) {
    Edge e(a,b);
    if (edges.find(e) != edges.end()) {
        // edge exists
        return;
    }
    edges.insert(e);
    Node na(a);
    Node nb(b);
    add_node(a);
    adj_map[na].insert(nb);
    adj_map[nb].insert(na);
    add_node(b);
    nodes.increment_degree(na);
    nodes.increment_degree(nb);
}

set<Node> Graph::neighbors(const Node& n) const {
    return adj_map.find(n)->second;//.second;
}

bool Graph::has_node(const Node &node) const { return nodes.exists(node); }

Graph* Graph::subgraph(const vector<Node> &nodes_container) const {
    Graph* G = new Graph();
    // add nodes
    for (Node node: nodes_container) {
        G->add_node(node);
    }
    // add edges
    for (Node node : G->get_nodes()) {
        // find neighbors in original graph
        for (Node neigh : this->neighbors(node)) {
            if (G->has_node(neigh)) {
                G->add_edge(node, neigh);
            }
        }
    }
    return G;
}

Graph* Graph::large_component() const {
    Components cc(*this);
    int max_component_id = cc.max_component_id();
    vector<Node> lcc_nodes;
    for (auto node : this->get_nodes()) {
        if (cc.node_component_id(node) == max_component_id) {
            lcc_nodes.push_back(node);
        }
    }
    Graph* G = this->subgraph(lcc_nodes);
    return G;
}

/*
set<Node> Graph::neighbors(const string& n) {
    return this->neighbors(Node(n));
}
*/

Components::Components(const Graph & G) : G(G), components_count(0) {
    //cout << "in Components" << endl;
    //cout << G.number_of_nodes() << endl;
    for (auto n : G.get_nodes()) {
        //cout << "\t" << n << endl;
        if (!is_labelled_node(n)) {
            ccR(n);
            components_count++;
        }
    }
}

bool Components::is_labelled_node(const Node &n) const {
    return component_ids.find(n) != component_ids.end();
}

void Components::ccR(const Node &n) {
    component_ids[n] = components_count;
    if (component_sizes.find(components_count) != component_sizes.end()) {
        component_sizes[components_count]++;
    } else {
        component_sizes[components_count] = 1;
    }

    // traverse from node using GFS
    stack<Node> dfs_stack;
    dfs_stack.push(n);
    while(!dfs_stack.empty()) {
        Node s = dfs_stack.top();
        dfs_stack.pop();
        for (auto neigh : G.neighbors(s)) {
            if (!is_labelled_node(neigh)) {
                component_ids[neigh] = components_count;
                component_sizes[components_count]++;
                dfs_stack.push(neigh);
            }
        }
    }
}

int Components::node_component_id(const Node &n) {
    if (is_labelled_node(n)) {
        return component_ids[n];
    } else {
        return -1;
    }
}

unsigned int Components::component_size(int component_id) {
    if (component_sizes.find(component_id) != component_sizes.end()) {
        return component_sizes[component_id];
    } else {
        return 0;
    }
}

int Components::get_components_count() {return components_count;}

int Components::max_component_id() {
    int max_component_size = 0;
    int max_component_id = -1;
    for (auto it: component_sizes) {
        if (it.second >= max_component_size) {
            max_component_id = it.first;
            max_component_size = component_sizes[it.first];
        }
    }
    return max_component_id;
}