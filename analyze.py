#!/usr/bin/env python

import sys

from argparse import ArgumentParser
from collections import Counter

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

def parse_vertex(raw_vertex, layout):
    if layout in ['polar', '2d']:
        x1, x2 = [float(rv) for rv in raw_vertex.split(',')]
        return (x1, x2)
    else:
        return raw_vertex

def main():
    parser = ArgumentParser()
    parser.add_argument('-f', help='file with graph (edges)')
    parser.add_argument('--clustering', help='compute clustering coefficients', action='store_true')
    parser.add_argument('--diameter', help='compute diameter', action='store_true')
    parser.add_argument('--pg', action='store_true', help='plot graph')
    parser.add_argument('--layout', help='graph layout (no, polar, 2d)', default='no')
    parser.add_argument('--pd', action='store_true', help='plot degree distribution')

    args = parser.parse_args()
    if args.f:
        graph_f = open(args.f)
    else:
        graph_f = sys.stdin

    print "read graph"
    graph = nx.Graph()
    for line in graph_f:
        if line.startswith('#'):
            continue
        items = line.rstrip().split()
        v1 = parse_vertex(items[0], args.layout)
        v2 = parse_vertex(items[1], args.layout)
        graph.add_edge(v1, v2)

    print "number of vertices: {}".format(graph.number_of_nodes())
    print "number of edges: {}".format(graph.number_of_edges())

    n = graph.number_of_nodes()

    # analyze degree distribution
    deg = [d[1] for d in graph.degree()]
    n_isolated_vertices = len([d for d in deg if d == 0])
    deg_dist = [i for i in Counter(deg).iteritems() if i[0] > 0]
    degs, counts = zip(*deg_dist)
    freqs = [float(c) / n for c in counts]

    # components
    print "n isolated vertices: {}".format(n_isolated_vertices)
    components = sorted(nx.connected_component_subgraphs(graph), key=len, reverse=True)
    print "n components: {}".format(len(components))
    print "n nontrivial components: {}".format(len(components) - n_isolated_vertices)
    print "largest compoment size: {}".format(len(components[0]))

    # clustering
    if args.clustering:
        print "Global clustering coefficient (transitivity): {0:.5f}".format(nx.transitivity(graph))
        print "Local clustering coefficient (average clustering): {0:.5f}".format(nx.average_clustering(graph))

    # diameter
    if args.diameter:
        Gcc=sorted(nx.connected_component_subgraphs(graph), key = len, reverse=True)
        G0=Gcc[0]
        print "Diameter of giant component: {}".format(nx.diameter(G0))

    # plots
    plots = [
        ['graph', args.pg],
        ['degrees', args.pd],
    ]
    active_plots = [pl for pl in plots if pl[1] is True]
    plot2idx = {pl[0]: pl_i+1 for pl_i, pl in enumerate(active_plots)}
    n_plots = len(plot2idx)

    # plot degree distribution
    if args.pd:
        plot_idx = plot2idx['degrees']
        ax_deg = plt.subplot(1, n_plots, plot_idx)
        plt.title('Degree distribution')
        ax_deg.loglog(degs, freqs, marker='o', linewidth=0)
        # fit in logarithmic scale
        # truncate 10% boundary degrees
        x = np.log(degs)
        y = np.log(freqs)
        xy = sorted(zip(x, y), key=lambda z: z[0])
        truncate_xy = int(0.1 * len(xy))
        xy = xy[truncate_xy:-truncate_xy]
        x, y = zip(*xy)

        fit = np.polyfit(x, y, 1)
        ax_deg.loglog(np.exp(x), np.exp(np.poly1d(fit)(x)), label="{0:.2f}x{1:+.2f}".format(*fit))
        ax_deg.grid(True)
        ax_deg.legend()

    # plot graph
    if args.pg:
        if args.layout == 'polar':
            plot_idx = plot2idx['graph']
            # vertices
            r, phi = zip(*graph.nodes())
            ax = plt.subplot(1, n_plots, plot_idx, projection='polar')
            ax.plot(phi, r, marker='o', linewidth=0)

            # edges
            for (v1, v2) in graph.edges():
                r1, phi1 = v1
                r2, phi2 = v2
                ax.plot((phi1, phi2), (r1, r2), color='g')
        elif args.layout == '2d':
            raise Exception('2d layout is not implemented')
        else:
            raise Exception('default layout is not implemented')
    if n_plots > 0:
        plt.show()

if __name__ == '__main__':
    main()
