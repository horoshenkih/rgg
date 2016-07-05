#!/usr/bin/env python
import argparse
from collections import defaultdict
from itertools import combinations

import numpy as np

def cosh_d(v1, v2):
    r1, phi1 = v1
    r2, phi2 = v2

    return np.cosh(r1) * np.cosh(r2) - np.sinh(r1) * np.sinh(r2) * np.cos(phi1 - phi2)

def find_embeddings(vertices, edges):
    "find (r, phi) for each vertex"
    n = len(vertices)
    R = 2 * np.log(n)
    # benchmark: phi=rand, r = 2log(n/k)
    degrees = defaultdict(int)
    for v1, v2 in edges:
        degrees[v1] += 1
        degrees[v2] += 1

    np.random.seed(0)
    return {v: (2*np.log(n / degrees[v]), np.random.uniform(0.0, 2*np.pi)) for v in vertices}

def evaluate_embeddings(embeddings, edges):
    "evaluate quality of embeddings compared with real elge set"
    n = len(embeddings.keys())
    R = 2 * np.log(n)
    coshR = np.cosh(R)

    report = defaultdict(int)
    for v1, v2 in combinations(embeddings.keys(), 2):
        is_true_edge = (v1, v2) in edges or (v2, v1) in edges
        is_predicted_edge = cosh_d(embeddings[v1], embeddings[v2]) <= coshR
        if is_true_edge:
            if is_predicted_edge:
                report['true_positive'] += 1
            else:
                report['false_negative'] += 1
        else:
            if is_predicted_edge:
                report['false_positive'] += 1
            else:
                report['true_negative'] += 1

    return report

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('graph_file')

    args = parser.parse_args()
    vertices = set()
    edges = set()
    with open(args.graph_file, 'r') as f:
        for line in f:
            v1, v2 = line.rstrip().split()
            vertices.update((v1, v2))
            edges.add((v1, v2))
    n = len(vertices)
    print "Number of vertices: {}".format(n)
    print "Number of edges: {}".format(len(edges))
    print "Number of non-edges: {}".format(n*(n-1)/2 - len(edges))

    embeddings = find_embeddings(vertices, edges)
    print embeddings
    report = evaluate_embeddings(embeddings, edges)

    print 'report'
    for k in sorted(report.keys()):
        print '{}: {}'.format(k, report[k])

if __name__ == '__main__':
    main()
