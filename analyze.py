#!/usr/bin/env python

import sys

from argparse import ArgumentParser
from collections import Counter

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from lib.graph import read_embeddings_from_file, read_communities_from_file

def deg2rad(x):
    # faster to write than to google it
    return 2 * np.pi * x / 360.

def main():
    parser = ArgumentParser()
    parser.add_argument('-f', help='file with graph, edge in each line')

    parser.add_argument('-c', help='file with communities, list of vertices in each line')
    parser.add_argument('--n-top-comm', help='number of top communities', type=int)

    parser.add_argument('-e', help='file with embeddings, triple (vertex, r, phi) in each line')
    parser.add_argument('--emb-skip-lines', help='skip first lines in embedding file', type=int, default=0)
    parser.add_argument('--deg', action='store_true', help='angle in degrees, not radians')
    parser.add_argument('--annotate', action='store_true', help='annotate vertices')

    parser.add_argument('--clustering', help='compute clustering coefficients', action='store_true')

    parser.add_argument('--diameter', help='compute diameter', action='store_true')

    parser.add_argument('--pg', action='store_true', help='plot graph')
    parser.add_argument('--layout', help='graph layout (no, polar, 2d)', default='polar')
    parser.add_argument('--no-edges', action='store_true', help='do not plot graph edges')
    parser.add_argument('--n-vertices', type=int, help='plot number of vertices closest to the origin')
    parser.add_argument('--core', action='store_true', help='plot only core vertices')
    parser.add_argument('--core_exponent', type=float, default=0.5)

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

    communities = None
    if args.c:
        print "read communities"
        communities = read_communities_from_file(args.c, args.n_top_comm)

    embeddings = None
    if args.e:
        print "read embeddings"
        embeddings = read_embeddings_from_file(args.e, args.emb_skip_lines)

    print "number of vertices: {}".format(graph.number_of_nodes())
    print "number of edges: {}".format(graph.number_of_edges())

    n = graph.number_of_nodes()

    # analyze degree distribution
    deg = graph.degree().values()# [d[1] for d in graph.degree()]
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
            elif args.core:
                n = graph.number_of_nodes()
                vert = [v for v in vert if graph.degree(v) >= n**args.core_exponent]
                print "Number of core vertices: {}".format(len(vert))
            pairs = [embeddings[v] for v in vert]
            r, phi = zip(*pairs)
            if args.deg:
                phi = map(deg2rad, phi)

            colors = None
            if communities is not None:
                colors=[]
                # find dominating community
                for v in vert:
                    comm_index, comm_size, comm_color = max(communities.get(v, [(0, 1, (0.9,0.9,0.9, 0.25))]), key=lambda x: x[1])
                    colors.append(comm_color)
            ax = plt.subplot(1, n_plots, plot_idx, projection='polar')
            # edges
            vert = set(vert)  # check inclustion in set is faster
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
            if colors is None:
                ax.plot(phi, r, marker='o', linewidth=0)
            else:
                for phi, r, color in zip (phi, r, colors):
                    ax.scatter(phi, r, color=color)
            if args.annotate:
                for v in vert:
                    r, phi = embeddings[v]
                    if args.deg:
                        phi = deg2rad(phi)
                    ax.annotate(
                        str(v), xy=(phi, r), xytext=(phi, r),
                        horizontalalignment='left', verticalalignment='top'
                    )

        elif args.layout == '2d':
            raise Exception('2d layout is not implemented')
        else:
            raise Exception('default layout is not implemented')
    if n_plots > 0:
        plt.show()

if __name__ == '__main__':
    main()
