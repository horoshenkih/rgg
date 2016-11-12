#include <iostream>
#include <string>
#include "boost/program_options.hpp"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char **argv) {
    string graph_file;
    string out_prefix;

    // parse options
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("graph-file", po::value<string>(&graph_file)->required(), "file with graph nodes")
            ("out-prefix", po::value<string>(&out_prefix)->required(), "prefix for output files")
            ("help", "produce help message");
    po::positional_options_description pos_desc;
    pos_desc.add("graph-file", 1);
    pos_desc.add("out-prefix", 2);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);

    if (vm.count("help")) {
        cout << desc << endl;
        return 1;
    }
    po::notify(vm);

    cout << "Hello, World!" << endl;
    return 0;
}