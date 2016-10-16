#!/usr/bin/env python
import argparse
from collections import defaultdict
from itertools import combinations
import time
import datetime

import numpy as np
import networkx as nx
from sklearn.neighbors import BallTree

from lib.graph import read_graph_from_file, read_embeddings_from_file, distance, make_edge

def evaluate_embeddings(embeddings, edges, cda=True, greedy_routing=False):
    "evaluate quality of embeddings compared with real elge set"
    report = []

    # get connected component
    true_graph = nx.Graph()
    true_graph.add_edges_from(edges)
    Gcc=sorted(nx.connected_component_subgraphs(true_graph), key=len, reverse=True)
    true_graph=Gcc[0]

    # use BallTree for efficient graph construction
    print "construct BallTree"
    vertices = list(true_graph.nodes())
    n = len(vertices)
    embeddings_array = np.array([embeddings[v] for v in vertices])
    bt = BallTree(embeddings_array, metric=distance)

    degrees = defaultdict(int)
    print "compute number of correct directed arcs"
    for v1, v2 in edges:
        degrees[v1] += 1
        degrees[v2] += 1

    # compute number of correct DIRECTED arcs assuming that degrees are known
    if cda:
        max_vertices = 6000
        all_correct_arcs = set()
        cda_vertices = vertices[:]
        if len(cda_vertices) > max_vertices:
            np.random.shuffle(cda_vertices)
            cda_vertices = cda_vertices[:max_vertices]
        for v_i, v in enumerate(cda_vertices):
            start = time.time()
            degree = degrees[v]
            dist, ind = bt.query(embeddings[v], k=degree+1)  # one of neighbors is vertex inself
            neigh = [vertices[i] for i in ind[0].tolist() if vertices[i] != v]
            for ne in neigh:
                if make_edge(v, ne) in edges:
                    all_correct_arcs.add((v, ne))
            finish = time.time()
            print "DEBUG: {} / {}, time={}s".format(v_i + 1, len(cda_vertices), datetime.timedelta(seconds=finish-start))
        report.append(['ratio of correct arcs for known degrees', float(len(all_correct_arcs)) / (2 * len(edges))])

    if greedy_routing:
        print "compute greedy routing efficiency"
        max_pairs = 10000
        random_pairs = set()
        if n * (n-1) / 2 <= max_pairs:
            random_pairs = set(combinations(vertices, 2))
        else:
            while(len(random_pairs) < max_pairs):
                v1 = np.random.choice(vertices)
                v2 = np.random.choice(vertices)
                if v1 != v2:
                    random_pairs.add((v1, v2))

        total_distribution = defaultdict(int)
        success_distribution = defaultdict(int)
        for i, pair in enumerate(random_pairs):
            src, dst = pair
            # best path
            best_path_length = nx.shortest_path_length(true_graph, source=src, target=dst)
            total_distribution[0] += 1
            total_distribution[best_path_length] += 1

            # greedy path
            curr_src = src
            path_length = 0
            seen = set()
            while curr_src != dst:
                seen.add(curr_src)
                # find neighbor closest to destination
                unseen_neighbors = filter(lambda x: x not in seen, true_graph.neighbors(curr_src))
                if not len(unseen_neighbors):
                    # greedy algorithm stuck in 'leaf'
                    path_length = np.nan
                    break

                def curr_distance(v):
                    return distance(embeddings[dst], embeddings[v])
                closest_neigh = min(unseen_neighbors, key=curr_distance)
                path_length += 1
                curr_src = closest_neigh

            if path_length == best_path_length:
                success_distribution[0] += 1
                success_distribution[best_path_length] += 1
        all_success = success_distribution[0]
        all_total = total_distribution[0]
        all_ratio = float(all_success) / all_total * 100
        print "Total: {} / {} ({:.2f} %)".format(all_success, all_total, all_ratio)
        for pl in sorted(set(total_distribution.keys()) | set(success_distribution.keys())):
            if pl == 0:
                continue
            total = total_distribution.get(pl, 0)
            success = success_distribution.get(pl, 0)
            ratio = float(success) / total * 100
            print "Path length = {}: {} / {} ({:.2f} %)".format(pl, success, total, ratio)

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

    return report

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('graph_file')
    parser.add_argument('embeddings_file')
    parser.add_argument('--deg', action='store_true', help='angle in degrees')
    parser.add_argument('--cda', action='store_true', help='correct directed arcs')
    parser.add_argument('--gr', action='store_true', help='evaluate greedy routing')
    args = parser.parse_args()

    vertices, edges = read_graph_from_file(args.graph_file)
    skip_lines=0
    if args.deg:
        skip_lines=2
    embeddings = read_embeddings_from_file(args.embeddings_file, skip_lines, degrees=args.deg)

    print "Evaluate embeddings"
    report = evaluate_embeddings(embeddings, edges, cda=args.cda, greedy_routing=args.gr)

    print 'report'
    for name, value in report:
        print '{}: {}'.format(name, value)

if __name__ == '__main__':
    main()
