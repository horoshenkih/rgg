#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iterator>

#include "lib/utils.h"
#include "lib/pair_generator.h"
#include "lib/embedding_model.h"
#include "lib/loss_function.h"
#include "lib/optimization.h"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;

PoincareModel fit(Graph *G) {
    // construct connected (TODO!) core
    const double core_exponent = 0.5;
    set<Node> core = G->core_nodes(core_exponent);
    cout << "Core size: " << core.size() << endl;
    cout << "Is core subgraph connected? " << G->subgraph(core)->is_connected() << endl;
    PairGenerator pair_generator(G);
    cout << "Prepare embedding model" << endl;
    PoincareModel embedding(G);
    LogLoss loss_function;
    SGD optimizer(0.1, 10, true); // learning rate 0.2

    cout << "Start embedding optimization" << endl;
    optimizer.optimize_embedding(&embedding, &loss_function, &pair_generator, G);
    return embedding;

    /*
    cout << "Pairs:" << endl;
    for (auto e : pair_generator.get_pairs()) {
        cout << e.to_string() << endl;
    }
    */


    /*
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
     */
}

void describe_graph(Graph &G) {
    cout << "Number of nodes: " << G.number_of_nodes() << endl;
    cout << "Number of edges: " << G.number_of_edges() << endl;

    cout << "Find components" << endl;
    Components cc(G);
    cout << "Number of components: " << cc.get_components_count() << endl;
    /*
    for (const auto &n: G.get_sorted_nodes()) {
        int cid = cc.node_component_id(n);
        cout << n << "\t" << cid << "\t" << cc.component_size(cid) << "\t" << G.get_node_description(n).get_degree() << endl;
    }
     */
    int max_component_id = cc.max_component_id();
    cout << "Max component: " << max_component_id << endl;
    cout << "Max component size: " << cc.component_size(max_component_id) << endl;
}

int main(int argc, char **argv) {
    std::srand(42);
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
    PoincareModel embedding = fit(subG);
    cout << "Save embedings to " << out_embeddings << endl;
    write_embedding_to_file(embedding, out_embeddings.c_str());
    /*
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
*/
    return 0;
}
