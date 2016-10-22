#!/usr/bin/python

import sys
sys.path.append('..')
from lib.graph import read_graph_from_file
import networkx as nx
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('graph')

    args = parser.parse_args()
    #vertices, edges = read_graph_from_file(args.graph)
    G = nx.read_edgelist(args.graph)
    n = G.number_of_nodes()
    print "nodes:", n
    print "edges:", G.number_of_edges()
    core_vertices = filter(lambda v: G.degree(v) >= n**0.5, G.nodes())
    print "core vertices:", len(core_vertices)
    core = G.subgraph(core_vertices)
    print "number of connected components:", nx.number_connected_components(core)

if __name__ == '__main__':
    main()
