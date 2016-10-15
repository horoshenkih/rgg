import networkx as nx
import numpy as np
from scipy.stats import rv_discrete, rv_continuous, uniform
from collections import Counter

class hyp_radius_density(rv_continuous):
    def __init__(self, alpha, R):
        self.alpha = alpha
        self.R = R
        super(hyp_radius_density, self).__init__(a=0., b=R)

    def _pdf(self, r):
        return self.alpha * np.sinh(self.alpha * r) / (np.cosh(self.alpha * self.R) - 1.)

def hyp_distance_cosh(x1, x2):
    r1, phi1 = x1
    r2, phi2 = x2
    return np.cosh(r1) * np.cosh(r2) - np.sinh(r1) * np.sinh(r2) * np.cos(phi1 - phi2)

def hyp_distance(x1, x2):
    return np.arccosh(hyp_distance_cosh(x1, x2))

# hyperbolic random graph
# vertex is a tuple (r, phi)
# TODO inherit nx.Graph()?
class HypRG:
    def __init__(self, n, alpha=1., C=0., seed=0):
        self.n = n
        self.alpha = alpha
        self.C = C
        np.random.seed(seed)

        # prepare
        self.R = 2 * np.log(self.n) + self.C
        assert self.R > 0
        self.angle_distribution = uniform(loc=0., scale=2 * np.pi).rvs
        self.radius_distribution = hyp_radius_density(self.alpha, self.R)().rvs  # frozen

        self.graph = nx.Graph()
        # generate vertices
        self._vertices = self.generate_vertices()
        for i, v in enumerate(self._vertices):
            self.graph.add_node(v)
            if i > 0:
                for old_v in self._vertices[:i-1]:
                    if hyp_distance(old_v, v) < self.R:
                        self.graph.add_edge(old_v, v)

    def vertices(self):
        return self._vertices

    def edges(self):
        return list(self.graph.edges())

    def generate_vertices(self):
        # angle is uniform in [0, 2 * pi)
        # radius has density alpha * sinh(alpha * r) / (cosh(alpha * R) - 1)
        r = self.radius_distribution(size=self.n)
        phi = self.angle_distribution(size=self.n)
        return zip(r, phi)

    def degrees(self):
        return [d[1] for d in self.graph.degree()]

    def components(self):
        return sorted(nx.connected_component_subgraphs(self.graph), key=len, reverse=True)

# t = number of vertices
# n = dimension (integer)
# cluster_attachment = 'preferential' / 'uniform'
# new_vertex = 'density' / 'uniform'
class RGG:
    def __init__(self, t=1, n=1, cluster_attachment='preferential', new_vertex = 'density', seed=None):
        assert t > 0
        assert n == 1  # TODO
        self.n = n
        self.cluster_attachment = cluster_attachment
        self.new_vertex = new_vertex
        np.random.seed(seed)

        self.graph = nx.Graph()
        # vertex is a pair (coords, cluster), coords is a tuple of floats
        v0 = (tuple([0.0 for i in range(self.n)]), 1)
        self.graph.add_node(v0)
        self.clusters = Counter({1: 1})

        # aux
        self.bbox = [-0.5, 0.5]  # TODO 1-dimensional
        for i in range(t-1):
            self.add_vertex()

    def add_vertex(self):
        # Step 1. Choose cluster
        clusters = sorted(self.clusters.keys())
        new_cluster = len(clusters) + 1
        if self.cluster_attachment == 'preferential':
            # probability of new cluster is 1 / (n_vertices + 1)
            norm = float(sum(self.clusters.values())+1)
            p_new = 1. / norm
            probas = [self.clusters[cl] / norm for cl in clusters]
        elif self.cluster_attachment == 'uniform':
            # probability of new cluster is 1 / (n_clusters + 1)
            norm = float(len(self.clusters.keys())+1)
            p_new = 1. / norm
            probas = [1. / norm for cl in clusters]
        # determine if new cluster appears
        if np.random.random() < p_new:
            self.clusters[new_cluster] = 1
            # generate new vertex
            #coord = self._gen_1d_vertex(self.bbox[0], self.bbox[1])
            coord = np.random.choice(self.bbox)  # select one boundary
            new_vertex = ((coord,), new_cluster)
        else:
            distr = rv_discrete(values=(clusters, probas))
            cl = distr.rvs(size=1)[0]
            self.clusters[cl] += 1
            # generate vertex from existing cluster
            cluster_vertices = np.array([v[0][0] for v in self.graph.nodes() if v[1] == cl])
            pivot = np.random.choice(cluster_vertices)
            coord = self._gen_1d_vertex(pivot - 0.5, pivot + 0.5)
            new_vertex = ((coord,), cl)
        old_vertices = self.graph.nodes()
        self.graph.add_node(new_vertex)
        for v in old_vertices:
            if abs(v[0][0] - new_vertex[0][0]) < 0.5:
                self.graph.add_edge(v, new_vertex)

    def degrees(self):
        return [d[1] for d in self.graph.degree()]

    def _gen_1d_vertex(self, left, right):
        coord = np.random.uniform(left, right)
        new_bounds = [coord - 0.5, coord + 0.5]
        self.bbox[0] = min(self.bbox[0], new_bounds[0])
        self.bbox[1] = max(self.bbox[1], new_bounds[1])
        return coord

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

def grad_distance(v1, v2):
    cd = cosh_d(v1, v2)
    if abs(cd - 1.) < 1e-15:
        return np.array((0.,0.,0.,0.))
    grad_cd = grad_cosh_d(v1, v2)
    return grad_cd / np.sqrt(cd - 1) / np.sqrt(cd + 1)

def make_edge(v1, v2):
    return tuple(sorted((v1, v2)))

def read_graph_from_file(filename):
    with open(filename, 'r') as f:
        edges = set()
        vertices = set()
        for line in f:
            v1, v2 = line.rstrip().split()
            vertices.update((v1, v2))
            edges.add(make_edge(v1, v2))
        return vertices, edges

def read_embeddings_from_file(filename, skip_lines=0):
    with open(filename, 'r') as f:
        embeddings = dict()
        for i_line, line in enumerate(f):
            if i_line + 1 <= skip_lines:
                continue
            v, r, phi = line.rstrip().split()
            r = float(r)
            phi = float(phi)
            embeddings[v] = (r, phi)
        return embeddings

def read_communities_from_file(filename, n_top_communities=None):
    with open(filename, 'r') as f:
        communities = dict()
        all_comm_size = []
        np.random.seed(0)  # community colors should be the same if possible
        for i_line, line in enumerate(f):
            comm_index = i_line + 1
            vertices = line.rstrip().split()
            comm_size = len(vertices)
            all_comm_size.append((comm_index, comm_size))
            R = np.random.uniform()
            G = np.random.uniform()
            B = np.random.uniform()
            comm_color = (R,G,B)  # random RGB
            comm_info = (comm_index, comm_size, comm_color)
            for v in vertices:
                if v in communities:
                    communities[v].append(comm_info)
                else:
                    communities[v] = [comm_info]

        if n_top_communities is not None:
            # find top communities
            all_comm_size.sort(key=lambda x: -x[1])
            top_comms = set([comm for comm, size in all_comm_size[:n_top_communities]])

            # select only top communitites for each vertex
            for v, comm_infos in communities.items():
                top_comm_infos = filter(lambda x: x[0] in top_comms, comm_infos)
                if len(top_comm_infos):
                    communities[v] = top_comm_infos
                else:
                    del communities[v]
        return communities
