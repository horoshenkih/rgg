#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <iterator>

#include "lib/utils.h"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;

void describe_graph(const Graph &G) {
    cout << "Number of nodes: " << G.number_of_nodes() << endl;
    cout << "Number of edges: " << G.number_of_edges() << endl;

    cout << "Find components" << endl;
    Components cc(G);
    cout << "Number of components: " << cc.get_components_count() << endl;
    /*
    for (const auto &n: G.get_nodes()) {
        int cid = cc.node_component_id(n);
        cout << n << "\t" << cid << "\t" << cc.component_size(cid) << endl;
    }
     */
    int max_component_id = cc.max_component_id();
    cout << "Max component: " << max_component_id << endl;
    cout << "Max component size: " << cc.component_size(max_component_id) << endl;
}

void fit(Graph *G) {
    cout << "in fit!" << endl;
    // construct connected (!) core
    const double core_exponent = 0.5;
    set<Node> core = G->core_nodes(core_exponent);
    cout << "Core size: " << core.size() << endl;
    cout << "Is core subgraph connected? " << G->subgraph(core)->is_connected() << endl;
    cout << "Core nodes:" << endl;
    for (auto cn : core) {
        cout << "\t" << cn << endl;
    }

    set<Node> fringe = G->fringe_nodes(core);
    cout << "Fringe nodes:" << endl;
    for (auto fn : fringe) {
        cout << "\t" << fn << endl;
    }

    std::set_union(core.begin(), core.end(), fringe.begin(), fringe.end(), std::inserter(core, core.end()));
    cout << "Core nodes 2:" << endl;
    for (auto cn : core) {
        cout << "\t" << cn << endl;
    }
    cout << "Fringe nodes 2:" << endl;
    for (auto fn : G->fringe_nodes(core)) {
        cout << "\t" << fn << endl;
    }
}

int main(int argc, char **argv) {
    string graph_file;
    string out_prefix;

    // parse options
    vector<string> command_line(argv+1, argv+argc);
    if (command_line.size() < 2) {
        cout << "Usage: ./hyprg GRAPH_FILE OUT_PREFIX" << endl;
        return 1;
    }
    graph_file = command_line[0];
    out_prefix = command_line[1];

    string out_embeddings = out_prefix + "-" + "embeddings.txt";
    ofstream out_embeddings_file;
    out_embeddings_file.open(out_embeddings);

    cout << "Read graph from: " << graph_file << endl;
    Graph G = read_graph_from_file(graph_file.c_str());
    cout << "Original graph:" << endl;
    describe_graph(G);
    cout << "==========" << endl << endl;

    Graph* subG = G.large_component();
    cout << "Large component:" << endl;
    describe_graph(*subG);
    cout << "==========" << endl << endl;

    cout << "Fit" << endl;
    fit(subG);
    cout << "Save embedings to " << out_embeddings << endl;
    int n = subG->number_of_nodes();
    for (auto v : subG->get_nodes()) {
        NodeDescription d = subG->get_node_description(v);
        Coordinates c = d.get_coordinates();
        unsigned int deg = d.get_degree();
        double r = 2. * log(double(n) / deg);
        c[0] = r;
        float phi = 2 * 3.14 * static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        c[1] = phi;
        d.set_coordinates(c);
        subG->set_node_description(v, d);
        d = subG->get_node_description(v);
        c = d.get_coordinates();
        out_embeddings_file << v << "\t" << c[0] << "\t" << c[1] << endl;
    }
    out_embeddings_file.close();

    return 0;
}
