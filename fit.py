#!/usr/bin/env python
import argparse
from collections import defaultdict
from itertools import combinations
import os
import datetime
import time

import random
import numpy as np

from scipy.optimize import minimize, check_grad

import networkx as nx

from lib.graph import make_edge, read_graph_from_file, cosh_d, distance, grad_cosh_d

class Margin:
    "difference between distance and R"
    def __init__(self, R):
        self.R = R
        self.coshR = np.cosh(R)
    def __call__(self, r1, phi1, r2, phi2):
        cd = cosh_d((r1, phi1), (r2, phi2))
        return np.arccosh(cd) - self.R

class GradMargin(Margin):
    "gradient of margin wrt r1, phi1, r2, phi2"
    def __call__(self, r1, phi1, r2, phi2):
        cd = cosh_d((r1, phi1), (r2, phi2))
        if abs(cd - 1.) < 1e-15:
            return np.array((0.,0.,0.,0.))
        grad_cd = grad_cosh_d((r1, phi1), (r2, phi2))
        return grad_cd / np.sqrt(cd - 1) / np.sqrt(cd + 1)

class Smooth:
    """approximation of step function of margin
    also, model of edge probability"""
    def __init__(self, beta=1., height=1.):
        self.height = height
        self.beta = beta
    def __call__(self, margin):
        return 1 / (1. + np.exp(margin * self.beta)) * self.height

class GradSmooth(Smooth):
    "gradient of step function approximation wrt margin"
    def __call__(self, margin):
        return -np.exp(margin * self.beta) / (1. + np.exp(margin * self.beta))**2 * self.height * self.beta

class Q:
    "loss function to minimize"
    def __init__(self, vertices, edges, nedges):
        self.vertices = vertices
        self.edges = edges
        self.nedges = nedges
        n = len(vertices)
        assert n > 1
        R = 2 * np.log(n)
        self.R = R
        self.coshR = np.cosh(R)
        self.non_edge_weight = float(len(self.edges)) / len(self.nedges) if len(self.nedges) else 1.

        self.margin = Margin(R)
        self.grad_margin = GradMargin(R)

        beta = 1.
        self.smooth = Smooth(beta=beta)
        self.grad_smooth = GradSmooth(beta=beta)

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

    def vertex_pair_grad(self, x, v1, v2, is_true_edge):
        assert len(x) % 2 == 0

        value = np.zeros(len(x))
        i1 = self.vertices.index(v1)
        i2 = self.vertices.index(v2)
        edge_ind = 1. if is_true_edge else 0.
        v = self.__grad_terms(x, i1, i2, edge_ind)
        value[2*i1]   = v[0]  # r1
        value[2*i1+1] = v[1]  # phi1
        value[2*i2]   = v[2]  # r2
        value[2*i2+1] = v[3]  # phi2
        return value

def find_embeddings(vertices, edges, mode,
    learning_rate=0.1, n_epoch=100,
    ratio_to_second=2., ratio_between_first=1., ratio_random=1.):
    "find (r, phi) for each vertex"
    vertices = list(vertices)
    n = len(vertices)
    R = 2 * np.log(n)

    print "mode: {}".format(mode)
    np.random.seed(0)
    degrees = defaultdict(int)
    print "count degrees"
    for v1, v2 in edges:
        degrees[v1] += 1
        degrees[v2] += 1
    if mode=='random':
        # phi=rand(0, 2pi), r = rand(0,R)
        return {v: (np.random.uniform(0.0, R), np.random.uniform(0.0, 2*np.pi)) for v in vertices}
    elif mode == 'degrees':
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
            #if (v1, v2) not in edges and (v2, v1) not in edges:
            e = make_edge(v1, v2)
            if e not in edges:
                all_nedges.add(e)

        if mode == 'fit_random':
            a = list(all_nedges)
            random.shuffle(a)
            nedges = set(a[:len(edges)])
        elif mode.startswith('fit_degrees'):
            K = float(ratio_to_second)      # ratio of nedges to second neighbour
            L = float(ratio_between_first)  # ratio of nedges between first neighbours
            M = float(ratio_random)         # ratio of random nedges
            #free_nedges = all_nedges.copy()

            G = nx.Graph()
            G.add_edges_from(edges)
            srt_vertices = sorted(degrees.keys(), key=lambda v: -degrees[v])
            shuf_vertices = srt_vertices[:]
            random.shuffle(shuf_vertices)
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
                    e = make_edge(v, sec_n)
                    if e not in nedges:
                        nedges.add(e)
                        n_vertex_nedges += 1

                # between first neighbours
                for j, pair in enumerate(combinations(first_neigh, 2)):
                    #print "j: {}".format(j)
                    if j+1 > degrees[v] * L:
                        continue
                    v1, v2 = pair
                    e = make_edge(v1, v2)
                    if e not in nedges:
                        nedges.add(e)

                # random edges
                max_n_random_vertices = int(degrees[v]*M)
                n_random_vertices = 0
                for rand_v in shuf_vertices:
                    if n_random_vertices >= max_n_random_vertices:
                        break
                    e = make_edge(v, rand_v)
                    if e not in nedges and e not in edges:
                        nedges.add(e)
                        n_random_vertices += 1
        else:
            nedges = all_nedges.copy()
        print "number of nedges={}".format(len(nedges))
        q = Q(vertices, edges, nedges)
        grad_q = GradQ(vertices, edges, nedges)
        if mode == 'fit_degrees_sgd':
            print "Learning rate: {}".format(learning_rate)
            print "Ratio to second: {}".format(ratio_to_second)
            print "Ratio between first: {}".format(ratio_between_first)
            print "Ratio random: {}".format(ratio_random)
            x = x0
            triples = [(v1, v2, True) for v1, v2 in edges] + [(v1, v2, False) for v1, v2 in nedges]
            random.shuffle(triples)
            for epoch in range(n_epoch):
                print "Epoch {} / {} ...".format(epoch+1, n_epoch)
                start = time.time()
                for v1, v2, is_true_edge in triples:
                    x -= grad_q.vertex_pair_grad(x, v1, v2, is_true_edge) * learning_rate
                finish = time.time()
                print "Elapsed time: {}s".format(datetime.timedelta(seconds=finish-start))
        else:
            print "Check gradient: ", check_grad(q, grad_q, x0)
            res = minimize(q, x0, method='BFGS', jac=grad_q)
            #print res
            x = res.x
        retval = {}
        for i in range(len(vertices)):
            r = x[2*i]
            phi = x[2*i+1]
            retval[vertices[i]] = (r, phi)

        return retval
    else:
        raise Exception('unknown mode')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('graph_file')
    parser.add_argument('embeddings_outfile')
    parser.add_argument('--mode', default='fit', help='random|degrees|fit|fit_random|fit_degrees|fit_degrees_sgd')
    parser.add_argument('--learning_rate', default=0.1, help='learning rate for fit_degrees_sgd', type=float)
    parser.add_argument('--n_epoch', default=100, help='number of training epoch for fit_degrees_sgd', type=int)
    parser.add_argument('--ratio_to_second', default=2., help='ratio of nedges to second neighbour', type=float)
    parser.add_argument('--ratio_between_first', default=1., help='ratio of nedges between first neighbours', type=float)
    parser.add_argument('--ratio_random', default=1., help='ratio of random nedges', type=float)

    args = parser.parse_args()
    vertices, edges = read_graph_from_file(args.graph_file)
    n = len(vertices)
    print "Number of vertices: {}".format(n)
    print "Number of edges: {}".format(len(edges))
    print "Number of non-edges: {}".format(n*(n-1)/2 - len(edges))

    print "Find embeddings"
    embeddings = find_embeddings(vertices, edges, mode=args.mode,
        learning_rate=args.learning_rate, n_epoch=args.n_epoch,
        ratio_to_second=args.ratio_to_second, ratio_between_first=args.ratio_between_first, ratio_random=args.ratio_random
    )

    with open(args.embeddings_outfile, 'w') as of:
        for v in embeddings.keys():
            r, phi = embeddings[v]
            of.write(' '.join(map(str, [v, r, phi]))+'\n')

if __name__ == '__main__':
    main()
