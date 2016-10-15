#!/usr/bin/env python
import argparse
from collections import defaultdict

import numpy as np
import networkx as nx
from sklearn.neighbors import BallTree

from lib.graph import read_graph_from_file, read_embeddings_from_file, distance, make_edge

def evaluate_embeddings(embeddings, edges):
    "evaluate quality of embeddings compared with real elge set"
    report = []

    true_graph = nx.Graph()
    true_graph.add_edges_from(edges)

    # use BallTree for efficient graph construction
    print "construct BallTree"
    vertices = list(true_graph.nodes())
    embeddings_array = np.array([embeddings[v] for v in vertices])
    bt = BallTree(embeddings_array, metric=distance)

    if False:
        # depends on R, bad for subgraphs -- not used
        n = len(vertices)
        R = 2 * np.log(n)
        coshR = np.cosh(R)

        predicted_edges = set()
        print "predict edges"
        for v in vertices:
            coords = embeddings[v]
            neigh_idx = bt.query_radius(coords, R)
            neigh = [vertices[i] for i in neigh_idx[0].tolist() if vertices[i] != v]
            predicted_edges.update([make_edge(v, ne) for ne in neigh])
        report.append(['total_predicted_edges', len(predicted_edges)])

        # contingency matrix
        print "compute contingency matrix"
        report.append(['true positive', len(edges & predicted_edges)])
        report.append(['false positive', len(predicted_edges - edges)])
        report.append(['false negative', len(edges - predicted_edges)])
        report.append(['true negative', n*(n-1)/2 - len(edges | predicted_edges)])

    degrees = defaultdict(int)
    print "compute number of correct directed arcs"
    for v1, v2 in edges:
        degrees[v1] += 1
        degrees[v2] += 1

    # compute number of correct DIRECTED arcs assuming that degrees are known
    # TODO the same for non-edges
    all_correct_arcs = set()
    for v in vertices:
        degree = degrees[v]
        dist, ind = bt.query(embeddings[v], k=degree+1)  # one of neighbors is vertex inself
        neigh = [vertices[i] for i in ind[0].tolist() if vertices[i] != v]
        for ne in neigh:
            if make_edge(v, ne) in edges:
                all_correct_arcs.add((v, ne))
    report.append(['ratio of correct arcs for known degrees', float(len(all_correct_arcs)) / (2 * len(edges))])

    return report

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('graph_file')
    parser.add_argument('embeddings_file')
    parser.add_argument('-p', '--plot', action='store_true', help='plot visualization')
    args = parser.parse_args()

    vertices, edges = read_graph_from_file(args.graph_file)
    embeddings = read_embeddings_from_file(args.embeddings_file)

    print "Evaluate embeddings"
    report = evaluate_embeddings(embeddings, edges)

    print 'report'
    for name, value in report:
        print '{}: {}'.format(name, value)

    if args.plot:
        import matplotlib.pyplot as plt
        # vertices
        r, phi = zip(*embeddings.values())
        ax = plt.subplot(111, projection='polar')
        ax.plot(phi, r, marker='o', linewidth=0)

        # edges
        # TODO
        '''
        for (v1, v2) in graph.edges():
            r1, phi1 = v1
            r2, phi2 = v2
            ax.plot((phi1, phi2), (r1, r2), color='g')
        '''
        plt.show()

if __name__ == '__main__':
    main()
