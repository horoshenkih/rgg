#!/usr/bin/env python
import argparse
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.optimize import minimize

def cosh_d(v1, v2):
    r1, phi1 = v1
    r2, phi2 = v2

    return np.cosh(r1) * np.cosh(r2) - np.sinh(r1) * np.sinh(r2) * np.cos(phi1 - phi2)

class Q:
    "loss function to minimize"
    def __init__(self, vertices, edges):
        self.vertices = vertices
        self.edges = edges
        n = len(vertices)
        self.coshR = np.cosh(2 * np.log(n))

    def __call__(self, x):
        "x = [r1, phi1, r2, phi2, ...] for vertex sequence v1, v2, ..."
        value = 0.
        assert len(x) % 2 == 0
        for (v1, v2) in combinations(self.vertices, 2):
            i1 = self.vertices.index(v1)
            i2 = self.vertices.index(v2)
            r1 = x[2*i1]
            phi1 = x[2*i1+1]
            r2 = x[2*i2]
            phi2 = x[2*i2+1]
            pred_edge = 1. if cosh_d((r1, phi1), (r2, phi2)) <= self.coshR else 0.
            if (v1, v2) in self.edges or (v2, v1) in self.edges:
                true_edge = 1.
            else:
                true_edge = 0.
            value += (pred_edge - true_edge)**2
        return value

def find_embeddings(vertices, edges, mode):
    "find (r, phi) for each vertex"
    vertices = list(vertices)
    n = len(vertices)
    R = 2 * np.log(n)

    np.random.seed(0)
    if mode=='random':
        # phi=rand(0, 2pi), r = rand(0,R)
        return {v: (np.random.uniform(0.0, R), np.random.uniform(0.0, 2*np.pi)) for v in vertices}
    elif mode == 'degrees':
        # phi=rand(0,2pi), r = 2log(n/k)
        degrees = defaultdict(int)
        for v1, v2 in edges:
            degrees[v1] += 1
            degrees[v2] += 1
        return {v: (2*np.log(n / degrees[v]), np.random.uniform(0.0, 2*np.pi)) for v in vertices}
    elif mode == 'fit':
        x0 = []
        for (r, phi) in zip([np.random.uniform(0.0, R) for v in vertices], [np.random.uniform(0.0, 2*np.pi) for v in vertices]):
            x0.append(r)
            x0.append(phi)
        x0 = np.array(x0)

        q = Q(vertices, edges)
        res = minimize(q, x0)
        retval = {}
        for i in range(len(vertices)):
            r = res.x[2*i]
            phi = res.x[2*i+1]
            retval[vertices[i]] = (r, phi)

        return retval
    else:
        raise Exception('unknown mode')

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
    parser.add_argument('--mode', default='fit', help='random|degrees|fit')

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

    embeddings = find_embeddings(vertices, edges, mode=args.mode)
    report = evaluate_embeddings(embeddings, edges)

    print 'report'
    for k in sorted(report.keys()):
        print '{}: {}'.format(k, report[k])

if __name__ == '__main__':
    main()
