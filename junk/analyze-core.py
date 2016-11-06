#!/usr/bin/python

import sys
sys.path.append('..')
from lib.graph import read_graph_from_file, fringe
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
    core_exponent = 0.5
    core_vertices = filter(lambda v: G.degree(v) >= n**core_exponent, G.nodes())
    print "core vertices:", len(core_vertices)
    core = G.subgraph(core_vertices)
    print "number of connected components in core:", nx.number_connected_components(core)

    # BFS-traversal
    fringe_fraction = 0.1
    max_fringe_size = int(n * fringe_fraction)
    core_vertices = set(core_vertices)
    for i in range(int(1/fringe_fraction)+1):
        fringe_vertices = set(sorted(fringe(G, core_vertices), key=lambda v: -G.degree(v))[:max_fringe_size])
        if not fringe_vertices:
            break
        print "{}: core={}, fringe={}".format(i+1, len(core_vertices), len(fringe_vertices))
        core_vertices |= fringe_vertices

if __name__ == '__main__':
    main()
