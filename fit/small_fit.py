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
from sklearn.neighbors import BallTree, DistanceMetric

def make_edge(v1, v2):
    return tuple(sorted((v1, v2)))

def cosh_d(v1, v2):
    r1, phi1 = list(v1)[:2]
    r2, phi2 = list(v2)[:2]

    return max(1., np.cosh(r1) * np.cosh(r2) - np.sinh(r1) * np.sinh(r2) * np.cos(phi1 - phi2))  # Python precision issues

def distance(v1, v2):
    return np.arccosh(cosh_d(v1, v2))

def grad_cosh_d(v1, v2):
    r1, phi1 = v1
    r2, phi2 = v2

    ddr1 = np.sinh(r1) * np.cosh(r2) - np.cosh(r1) * np.sinh(r2) * np.cos(phi1 - phi2)
    ddr2 = np.cosh(r1) * np.sinh(r2) - np.sinh(r1) * np.cosh(r2) * np.cos(phi1 - phi2)
    ddphi1 = np.sinh(r1) * np.sinh(r2) * np.sin(phi1 - phi2)
    ddphi2 = -ddphi1

    return np.array((ddr1, ddphi1, ddr2, ddphi2))

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
            #if (v1, v2) not in edges and (v2, v1) not in edges:
            e = make_edge(v1, v2)
            if e not in edges:
                all_nedges.add(e)

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
    report = []

    n = len(embeddings.keys())
    R = 2 * np.log(n)
    coshR = np.cosh(R)

    #true_graph = nx.Graph()
    #true_graph.add_edges_from(edges)

    # use BallTree for efficient graph construction
    vertices = sorted(embeddings.keys())
    embeddings_array = np.array([embeddings[v] for v in vertices])
    bt = BallTree(embeddings_array, metric=distance)

    predicted_edges = set()
    for v in embeddings.keys():
        coords = embeddings[v]
        neigh_idx = bt.query_radius(coords, R)
        neigh = [vertices[i] for i in neigh_idx[0].tolist() if vertices[i] != v]
        predicted_edges.update([make_edge(v, ne) for ne in neigh])
    report.append(['total_predicted_edges', len(predicted_edges)])

    # contingency matrix
    report.append(['true positive', len(edges & predicted_edges)])
    report.append(['false positive', len(predicted_edges - edges)])
    report.append(['false negative', len(edges - predicted_edges)])
    report.append(['true negative', n*(n-1)/2 - len(edges | predicted_edges)])

    degrees = defaultdict(int)
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
    parser.add_argument('--mode', default='fit', help='random|degrees|fit|fit_random|fit_degrees')
    parser.add_argument('-p', '--plot', action='store_true', help='plot visualization')

    args = parser.parse_args()
    vertices = set()
    edges = set()
    with open(args.graph_file, 'r') as f:
        for line in f:
            v1, v2 = line.rstrip().split()
            vertices.update((v1, v2))
            edges.add(make_edge(v1, v2))
    n = len(vertices)
    print "Number of vertices: {}".format(n)
    print "Number of edges: {}".format(len(edges))
    print "Number of non-edges: {}".format(n*(n-1)/2 - len(edges))

    print "Find embeddings"
    embeddings = find_embeddings(vertices, edges, mode=args.mode)
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
        '''
        for (v1, v2) in graph.edges():
            r1, phi1 = v1
            r2, phi2 = v2
            ax.plot((phi1, phi2), (r1, r2), color='g')
        '''
        plt.show()

if __name__ == '__main__':
    main()
