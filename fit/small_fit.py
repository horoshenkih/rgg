#!/usr/bin/env python
import argparse
from collections import defaultdict
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, check_grad

def cosh_d(v1, v2):
    r1, phi1 = v1
    r2, phi2 = v2

    return np.cosh(r1) * np.cosh(r2) - np.sinh(r1) * np.sinh(r2) * np.cos(phi1 - phi2)

def grad_cosh_d(v1, v2):
    r1, phi1 = v1
    r2, phi2 = v2

    ddr1 = np.sinh(r1) * np.cosh(r2) - np.cosh(r1) * np.sinh(r2) * np.cos(phi1 - phi2)
    ddr2 = np.cosh(r1) * np.sinh(r2) - np.sinh(r1) * np.cosh(r2) * np.cos(phi1 - phi2)
    ddphi1 = np.sinh(r1) * np.sinh(r2) * np.sin(phi1 - phi2)
    ddphi2 = -ddphi1

    return np.array((ddr1, ddphi1, ddr2, ddphi2))

class Margin:
    "difference between log cosh(distance) and log cosh(R)"
    def __init__(self, R):
        self.coshR = np.cosh(R)
    def __call__(self, r1, phi1, r2, phi2):
        cd = cosh_d((r1, phi1), (r2, phi2))
        return np.log(cd / self.coshR)

class GradMargin(Margin):
    "gradient of margin wrt r1, phi1, r2, phi2"
    def __call__(self, r1, phi1, r2, phi2):
        cd = cosh_d((r1, phi1), (r2, phi2))
        grad_cd = grad_cosh_d((r1, phi1), (r2, phi2))
        return grad_cd / cd

class Smooth:
    "approximation of step function of margin"
    def __call__(self, margin):
        return 1 / (1. + np.exp(margin))

class GradSmooth(Smooth):
    "gradient of step function approximation wrt margin"
    def __call__(self, margin):
        return -np.exp(margin) / (1. + np.exp(margin))**2

class Q:
    "loss function to minimize"
    def __init__(self, vertices, edges):
        self.vertices = vertices
        self.edges = edges
        n = len(vertices)
        assert n > 1
        R = 2 * np.log(n)
        self.coshR = np.cosh(R)
        self.non_edge_weight = 2. * len(self.edges) / n / (n-1)

        self.margin = Margin(R)
        self.grad_margin = GradMargin(R)
        self.smooth = Smooth()
        self.grad_smooth = GradSmooth()

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

            z = self.margin(r1, phi1, r2, phi2)
            pred_edge = self.smooth(z)
            #pred_edge = 1. if cosh_d((r1, phi1), (r2, phi2)) <= self.coshR else 0.
            if (v1, v2) in self.edges or (v2, v1) in self.edges:
                true_edge = 1.
            else:
                true_edge = 0.
            w = 1. if true_edge else self.non_edge_weight
            value += (pred_edge - true_edge)**2 * w
        return value

class GradQ(Q):
    def __call__(self, x):
        assert len(x) % 2 == 0

        value = np.zeros(len(x))
        for (v1, v2) in combinations(self.vertices, 2):
            i1 = self.vertices.index(v1)
            i2 = self.vertices.index(v2)
            r1 = x[2*i1]
            phi1 = x[2*i1+1]
            r2 = x[2*i2]
            phi2 = x[2*i2+1]

            z = self.margin(r1, phi1, r2, phi2)
            smooth_der = self.grad_smooth(z)
            margin_der = self.grad_margin(r1, phi1, r2, phi2)
            true_edge = 1. if ((v1, v2) in self.edges or (v2, v1) in self.edges) else 0.
            w = 1. if true_edge else self.non_edge_weight
            disc = 2 * (self.smooth(z) - true_edge) * w
            value[2*i1]   += disc * smooth_der * margin_der[0]  # r1
            value[2*i1+1] += disc * smooth_der * margin_der[1]  # phi1
            value[2*i2]   += disc * smooth_der * margin_der[2]  # r2
            value[2*i2+1] += disc * smooth_der * margin_der[3]  # phi2

        return value

def find_embeddings(vertices, edges, mode):
    "find (r, phi) for each vertex"
    vertices = list(vertices)
    n = len(vertices)
    R = 2 * np.log(n)

    np.random.seed(0)
    degrees = defaultdict(int)
    for v1, v2 in edges:
        degrees[v1] += 1
        degrees[v2] += 1
    if mode=='random':
        # phi=rand(0, 2pi), r = rand(0,R)
        return {v: (np.random.uniform(0.0, R), np.random.uniform(0.0, 2*np.pi)) for v in vertices}
    elif mode == 'degrees':
        # phi=rand(0,2pi), r = 2log(n/k)
        return {v: (2*np.log(n / degrees[v]), np.random.uniform(0.0, 2*np.pi)) for v in vertices}
    elif mode == 'fit':
        x0 = []
        for (r, phi) in zip([2*np.log(n / degrees[v]) for v in vertices], [np.random.uniform(0.0, 2*np.pi) for v in vertices]):
            x0.append(r)
            x0.append(phi)
        x0 = np.array(x0)

        q = Q(vertices, edges)
        grad_q = GradQ(vertices, edges)
        print "Check gradient: ", check_grad(q, grad_q, x0)
        #res = minimize(q, x0, method='Nelder-Mead')
        #res = minimize(q, x0, method='Newton-CG', jac=grad_q)
        res = minimize(q, x0, method='BFGS', jac=grad_q)
        #print res
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
    parser.add_argument('-p', '--plot', action='store_true', help='plot visualization')

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
    if args.plot:
        # vertices
        r, phi = zip(*embeddings.values())
        ax = plt.subplot(111, projection='polar')
        ax.plot(phi, r, marker='o', linewidth=0)

        # edges
        '''
        for (v1, v2) in graph.edges():
            r1, phi1 = v1
            r2, phi2 = v2
            ax.plot((phi1, phi2), (r1, r2), color='g')
        '''
        plt.show()

if __name__ == '__main__':
    main()
