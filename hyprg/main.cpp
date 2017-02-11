#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iterator>

#include "lib/commandline.h"
#include "lib/utils.h"
#include "lib/pair_generator.h"
#include "lib/embedding_model.h"
#include "lib/loss_function.h"
#include "lib/optimization.h"

using std::cout;
using std::endl;
using std::string;
using std::ofstream;

void fit(Graph *G, PoincareModel* embedding,
         double learning_rate, int n_epoch,
         double ratio_to_second,
         double ratio_between_first,
         double ratio_first_second,
         double ratio_between_second,
         double ratio_random
) {
    // construct connected (TODO!) core
    const double core_exponent = 0.5;
    set<Node> core = G->core_nodes(core_exponent);
    cout << "Core size: " << core.size() << endl;
    cout << "Is core subgraph connected? " << G->subgraph(core)->is_connected() << endl;
    PairGenerator pair_generator(G, ratio_to_second, ratio_between_first, ratio_first_second, ratio_between_second, ratio_random);
    cout << "Prepare embedding model" << endl;
    //PoincareModel embedding(G);
    LogLoss loss_function;
    SGD optimizer(learning_rate, n_epoch, true);

    cout << "Start embedding optimization" << endl;
    optimizer.optimize_embedding(embedding, &loss_function, &pair_generator, G);
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
    CommandLine cl;

    cl.add_string_argument("g", "(required) file with graph edges");
    cl.add_string_argument("o", "(required) prefix of outfile");
    cl.add_string_argument("e", "file with embeddings");
    cl.add_double_argument("learning-rate", 0.2, "learning rate");
    cl.add_double_argument("ratio-to-second", 2., "ratio of nedges to second neighbours");
    cl.add_double_argument("ratio-between-first", 2., "ratio of nedges between first neighbours");
    cl.add_double_argument("ratio-first-second", 1., "ratio of nedges between first and second neighbours");
    cl.add_double_argument("ratio-between-second", 1., "ratio of nedges between second neighbours");
    cl.add_double_argument("ratio-random", 1., "ratio of nedges to random vertices");
    cl.add_int_argument("n-epoch", 20, "number of training epoch");
    cl.add_int_argument("seed", 42, "random seed");

    cl.parse_arguments(argc, argv);
    if (!cl.has_string_argument("g") or !cl.has_string_argument("o")) {
        cout << "Some options not provided, run program with -h or --help" << endl;
        return 1;
    }

    std::srand(cl.get_int_argument("seed"));
    string graph_file = cl.get_string_argument("g");
    string out_prefix = cl.get_string_argument("o");

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
    PoincareModel embedding(subG);
    if (cl.has_string_argument("e")) {
        cout << "Read initial ebmeddings from" << cl.get_string_argument("e") << endl;
        read_poincare_embedding_from_file(&embedding, cl.get_string_argument("e").c_str());
    }

    double learning_rate = cl.get_double_argument("learning-rate");
    int n_epoch = cl.get_int_argument("n-epoch");
    double ratio_to_second = cl.get_double_argument("ratio-to-second");
    double ratio_between_first = cl.get_double_argument("ratio-between-first");
    double ratio_first_second = cl.get_double_argument("ratio-first-second");
    double ratio_between_second = cl.get_double_argument("ratio-between-second");
    double ratio_random = cl.get_double_argument("ratio-random");

    fit(subG, &embedding, learning_rate, n_epoch, ratio_to_second, ratio_between_first, ratio_first_second, ratio_between_second, ratio_random);
    cout << "Save embedings to " << out_embeddings << endl;
    write_embedding_to_file(embedding, out_embeddings.c_str());

    return 0;
}
