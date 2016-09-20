from abc import ABCMeta, abstractmethod
from itertools import combinations
from collections import defaultdict
import random
import networkx as nx

from graph import make_edge

class PairGenerator:
    __metaclass__ = ABCMeta
    @abstractmethod
    def __init__(self, edges, **kwargs):
        pass

    @abstractmethod
    def __call__(self, **kwargs):
        pass

    def _count_degrees(self, edges):
        degrees = defaultdict(int)
        for v1, v2 in edges:
            degrees[v1] += 1
            degrees[v2] += 1
        return degrees

class BinaryPairGenerator(PairGenerator):
    def __init__(self, graph,
        ratio_to_second=2., ratio_between_first=1., ratio_random=1.,
        seed=0):
        random.seed(seed)
        self.graph = graph
        degrees = self._count_degrees(graph.edges())
        srt_vertices = sorted(degrees.keys(), key=lambda v: -degrees[v])
        shuf_vertices = srt_vertices[:]
        random.shuffle(shuf_vertices)

        edge_set = set(graph.edges())
        nedge_set = set()
        for v in srt_vertices:
            # get first neighbours
            first_neigh = set(self.graph.neighbors(v))
            # get second neighbours
            second_neigh = set()
            for neigh in first_neigh:
                second_neigh.update(self.graph.neighbors(neigh))
            second_neigh.remove(v)
            second_neigh = list(second_neigh)
            random.shuffle(second_neigh)

            n_second_vertex_nedges = 0
            # from v to second neighbours
            for sec_n in second_neigh:
                if n_second_vertex_nedges > degrees[v] * ratio_to_second:
                    break
                e = make_edge(v, sec_n)
                if e not in edge_set and e not in nedge_set:
                    nedge_set.add(e)
                    n_second_vertex_nedges += 1

            # between first neighbours
            n_between_first_nedges = 0
            for pair in combinations(first_neigh, 2):
                if n_between_first_nedges > degrees[v] * ratio_between_first:
                    break
                v1, v2 = pair
                e = make_edge(v1, v2)
                if e not in edge_set and e not in nedge_set:
                    nedge_set.add(e)
                    n_between_first_nedges += 1

            # random edges
            max_n_random_vertices = int(degrees[v]*ratio_random)
            n_random_vertices = 0
            for rand_v in shuf_vertices:
                if rand_v == v:
                    continue
                if n_random_vertices >= max_n_random_vertices:
                    break
                e = make_edge(v, rand_v)
                if e not in edge_set and e not in nedge_set:
                    nedge_set.add(e)
                    n_random_vertices += 1

        self.pairs = [(e, True) for e in edge_set] + [(ne, False) for ne in nedge_set]
        random.shuffle(self.pairs)

    def __call__(self):
        for e in self.pairs:
            yield e

if __name__ == '__main__':
    edges = [
        (0,1), (0,2),
        (1, 11), (1, 12), (2, 21), (2, 22),
    ]
    G = nx.Graph()
    G.add_edges_from(edges)
    pair_gen = BinaryPairGenerator(G)
    for e in pair_gen():
        print e
