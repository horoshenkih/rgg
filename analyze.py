#!/usr/bin/env python

import sys

from argparse import ArgumentParser
from collections import Counter

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from lib.graph import read_embeddings_from_file

def main():
    parser = ArgumentParser()
    parser.add_argument('-f', help='file with graph (edges)')

    parser.add_argument('-e', help='file with embedding (vertex, r, phi)')
    parser.add_argument('--skip-lines', help='skip first lines in embedding file', type=int, default=0)
    parser.add_argument('--deg', action='store_true', help='angle in degrees, not radians')

    parser.add_argument('--clustering', help='compute clustering coefficients', action='store_true')

    parser.add_argument('--diameter', help='compute diameter', action='store_true')

    parser.add_argument('--pg', action='store_true', help='plot graph')
    parser.add_argument('--layout', help='graph layout (no, polar, 2d)', default='polar')
    parser.add_argument('--no-edges', action='store_true', help='do not plot graph edges')
    parser.add_argument('--n-vertices', type=int, help='plot number of vertices closest to the origin')

    parser.add_argument('--pd', action='store_true', help='plot degree distribution')

    args = parser.parse_args()
    if args.f:
        graph_f = open(args.f)
    else:
        graph_f = sys.stdin

    print "read graph"
    graph = nx.Graph()
    for i_line, line in enumerate(graph_f):
        if line.startswith('#'):
            continue
        v1, v2 = line.rstrip().split()
        graph.add_edge(v1, v2)

    embeddings = None
    if args.e:
        print "read embeddings"
        embeddings = read_embeddings_from_file(args.e, args.skip_lines)

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
            if embeddings is None:
                raise Exception("Embedding is not provided for polar layout")
            plot_idx = plot2idx['graph']

            vert = embeddings.keys()
            if args.n_vertices:
                vert.sort(key=lambda x: embeddings[x][0])
                vert = vert[:args.n_vertices]
            vert = set(vert)
            pairs = [embeddings[v] for v in vert]
            r, phi = zip(*pairs)
            def deg2rad(x):
                return 2 * np.pi * x / 360.
            if args.deg:
                phi = map(deg2rad, phi)

            ax = plt.subplot(1, n_plots, plot_idx, projection='polar')
            # edges
            if not args.no_edges:
                for (v1, v2) in graph.edges():
                    if v1 not in vert or v2 not in vert:
                        continue
                    r1, phi1 = embeddings[v1]
                    r2, phi2 = embeddings[v2]
                    if args.deg:
                        phi1 = deg2rad(phi1)
                        phi2 = deg2rad(phi2)
                    ax.plot((phi1, phi2), (r1, r2), color='g')
            # vertices
            ax.plot(phi, r, marker='o', linewidth=0)

        elif args.layout == '2d':
            raise Exception('2d layout is not implemented')
        else:
            raise Exception('default layout is not implemented')
    if n_plots > 0:
        plt.show()

if __name__ == '__main__':
    main()
