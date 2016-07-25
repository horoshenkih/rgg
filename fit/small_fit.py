#!/usr/bin/env python
import argparse
from collections import defaultdict
from itertools import combinations
import os

import random
import numpy as np
from scipy.optimize import minimize, check_grad
import networkx as nx
from sklearn.metrics import roc_auc_score

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
    def __init__(self, scale=1., height=1.):
        self.height = height
        self.scale = scale
    def __call__(self, margin):
        return 1 / (1. + np.exp(margin * self.scale)) * self.height

class GradSmooth(Smooth):
    "gradient of step function approximation wrt margin"
    def __call__(self, margin):
        return -np.exp(margin * self.scale) / (1. + np.exp(margin * self.scale))**2 * self.height * self.scale

class Q:
    "loss function to minimize"
    def __init__(self, vertices, edges, nedges):
        self.vertices = vertices
        self.edges = edges
        self.nedges = nedges
        n = len(vertices)
        assert n > 1
        R = 2 * np.log(n)
        self.coshR = np.cosh(R)
        #self.non_edge_weight = 2. * len(self.edges) / n / (n-1)
        self.non_edge_weight = float(len(self.edges)) / len(self.nedges) if len(self.nedges) else 1.

        self.margin = Margin(R)
        self.grad_margin = GradMargin(R)

        smooth_scale = 1.
        self.smooth = Smooth(scale=smooth_scale)
        self.grad_smooth = GradSmooth(scale=smooth_scale)

    def __value_term(self, x, v1, v2, true_edge):
        i1 = self.vertices.index(v1)
        i2 = self.vertices.index(v2)
        r1 = x[2*i1]
        phi1 = x[2*i1+1]
        r2 = x[2*i2]
        phi2 = x[2*i2+1]

        z = self.margin(r1, phi1, r2, phi2)
        pred_edge = self.smooth(z)
        w = 1. if true_edge else self.non_edge_weight
        return (pred_edge - true_edge)**2 * w

    def __call__(self, x):
        "x = [r1, phi1, r2, phi2, ...] for vertex sequence v1, v2, ..."
        value = 0.
        assert len(x) % 2 == 0
        for (v1, v2) in self.edges:
            value += self.__value_term(x, v1, v2, 1.)
        for (v1, v2) in self.nedges:
            value += self.__value_term(x, v1, v2, 0.)
        return value

class GradQ(Q):
    def __grad_terms(self, x, i1, i2, true_edge):
        r1 = x[2*i1]
        phi1 = x[2*i1+1]
        r2 = x[2*i2]
        phi2 = x[2*i2+1]

        z = self.margin(r1, phi1, r2, phi2)
        smooth_der = self.grad_smooth(z)
        margin_der = self.grad_margin(r1, phi1, r2, phi2)

        v1 = self.vertices[i1]
        v2 = self.vertices[i2]
        w = 1. if true_edge else self.non_edge_weight
        disc = 2 * (self.smooth(z) - true_edge) * w

        return disc * smooth_der * margin_der

    def __call__(self, x):
        assert len(x) % 2 == 0

        value = np.zeros(len(x))
        for (v1, v2) in self.edges:
            i1 = self.vertices.index(v1)
            i2 = self.vertices.index(v2)
            v = self.__grad_terms(x, i1, i2, 1.)
            value[2*i1]   += v[0]  # r1
            value[2*i1+1] += v[1]  # phi1
            value[2*i2]   += v[2]  # r2
            value[2*i2+1] += v[3]  # phi2
        for (v1, v2) in self.nedges:
            i1 = self.vertices.index(v1)
            i2 = self.vertices.index(v2)
            v = self.__grad_terms(x, i1, i2, 0.)
            value[2*i1]   += v[0]  # r1
            value[2*i1+1] += v[1]  # phi1
            value[2*i2]   += v[2]  # r2
            value[2*i2+1] += v[3]  # phi2
        return value

def find_embeddings(vertices, edges, mode):
    "find (r, phi) for each vertex"
    vertices = list(vertices)
    n = len(vertices)
    R = 2 * np.log(n)

    np.random.seed(0)
    degrees = defaultdict(int)
    print "count degrees"
    for v1, v2 in edges:
        degrees[v1] += 1
        degrees[v2] += 1
    if mode=='random':
        print "mode: random"
        # phi=rand(0, 2pi), r = rand(0,R)
        return {v: (np.random.uniform(0.0, R), np.random.uniform(0.0, 2*np.pi)) for v in vertices}
    elif mode == 'degrees':
        print "mode: degrees"
        # phi=rand(0,2pi), r = 2log(n/k)
        return {v: (2*np.log(n / degrees[v]), np.random.uniform(0.0, 2*np.pi)) for v in vertices}
    elif mode.startswith('fit'):
        x0 = []
        for (r, phi) in zip([2*np.log(n / degrees[v]) for v in vertices], [np.random.uniform(0.0, 2*np.pi) for v in vertices]):
            x0.append(r)
            x0.append(phi)
        x0 = np.array(x0)

        nedges = set()
        all_nedges = set()
        for (v1, v2) in combinations(vertices, 2):
            if (v1, v2) not in edges and (v2, v1) not in edges:
                all_nedges.add((v1, v2))

        if mode == 'fit_random':
            a = list(all_nedges)
            random.shuffle(a)
            nedges = set(a[:len(edges)])
        elif mode == 'fit_degrees':
            print "mode: fit_degrees"
            K = 2.  # ratio of nedges to second neighbour
            L = 1.  # ratio of nedges between first neighbours
            M = 1.  # ratio of random nedges
            #free_nedges = all_nedges.copy()

            G = nx.Graph()
            G.add_edges_from(edges)
            srt_vertices = sorted(degrees.keys(), key=lambda v: -degrees[v])
            for v in srt_vertices:
                # get first neighbours
                first_neigh = set(G.neighbors(v))
                # get second neighbours
                second_neigh = set()
                for neigh in first_neigh:
                    second_neigh.update(G.neighbors(neigh))
                second_neigh.remove(v)

                n_vertex_nedges = 0
                # from v to second neighbours
                for i, sec_n in enumerate(second_neigh):
                    #print "i: {}".format(i)
                    if i+1 > degrees[v] * K:
                        continue
                    #if (v, sec_n) in free_nedges or (sec_n, v) in free_nedges:
                    if (v, sec_n) not in nedges and (sec_n, v) not in nedges:
                        nedges.add((v, sec_n))
                        n_vertex_nedges += 1
                        #try:
                        #    free_nedges.remove((v, sec_n))
                        #except KeyError:
                        #    free_nedges.remove((sec_n, v))

                # between first neighbours
                for j, pair in enumerate(combinations(first_neigh, 2)):
                    #print "j: {}".format(j)
                    if j+1 > degrees[v] * L:
                        continue
                    v1, v2 = pair
                    #if (v1, v2) in free_nedges or (v2, v1) in free_nedges:
                    if (v1, v2) not in nedges and (v2, v1) not in nedges:
                        nedges.add((v1, v2))
                        #try:
                        #    free_nedges.remove((v1, v2))
                        #except KeyError:
                        #    free_nedges.remove((v2, v1))

                # random edges
                #n_random_vertices = int(degrees[v]*M)
                #n_free_vertex_nedges = (n-1) - n_vertex_nedges
                #a = list(free_nedges)
                #random.shuffle(a)
                #random_to_update = set(a[:int(degrees[v]*M)])
                #nedges.update(random_to_update)
                #free_nedges.difference_update(random_to_update)
            print "fit_degrees: number of nedges={}".format(len(nedges))
        else:
            nedges = all_nedges.copy()
        q = Q(vertices, edges, nedges)
        grad_q = GradQ(vertices, edges, nedges)
        print "Check gradient: ", check_grad(q, grad_q, x0)
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
    # contingency matrix
    if os.environ['DEBUG']:
        print "DEBUG: evaluate contingency matrix"
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

    # vertexwise ROC-AUC
    if os.environ['DEBUG']:
        print "DEBUG: evaluate vertexwise AUC"
    all_vertex_aucs = []
    for i_v1, v1 in enumerate(embeddings.keys()):
        true = []
        predicted = []
        for v2 in embeddings.keys():
            is_true_edge = (v1, v2) in edges or (v2, v1) in edges
            is_predicted_edge = cosh_d(embeddings[v1], embeddings[v2])
            true.append(is_true_edge)
            predicted.append(-is_predicted_edge)
        auc = roc_auc_score(true, predicted)
        all_vertex_aucs.append(auc)
    report['vertexwise_auc'] = np.mean(all_vertex_aucs)
    return report

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('graph_file')
    parser.add_argument('--mode', default='fit', help='random|degrees|fit|fit_random|fit_degrees')
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

    print "Find embeddings"
    embeddings = find_embeddings(vertices, edges, mode=args.mode)
    print "Evaluate embeddings"
    report = evaluate_embeddings(embeddings, edges)

    print 'report'
    for k in sorted(report.keys()):
        print '{}: {}'.format(k, report[k])
    if args.plot:
        import matplotlib.pyplot as plt
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
