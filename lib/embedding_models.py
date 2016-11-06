from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import sparse

from graph import distance, grad_distance, fringe

class EmbeddingModel:
    __metaclass__ = ABCMeta
    @abstractmethod
    def __init__(self, graph, **kwargs):
        pass

    @abstractmethod
    def get_state_vector(self):
        """Get vector of model state.
        It may include coordinates of embeddings, radius (radii), etc.
        """
        pass

    @abstractmethod
    def set_state_vector(self, x):
        """Set model state from vector"""
        pass

    @abstractmethod
    def get_distance_info(self, edges_batch):
        """Get distances, gradients and other information for given set of edges"""
        pass

class PoincareModel(EmbeddingModel):
    def __init__(self, graph, init_embedding=None, radius=None, fit_radius=False):
        self.graph = graph
        n = self.graph.number_of_nodes()
        self.fit_radius = fit_radius
        if self.fit_radius:
            self.__expected_vector_size = 2 * n + 1
        else:
            self.__expected_vector_size = 2 * n

        # initialize embedding
        if radius is None:
            radius = 2 * np.log(n)

        self.embedding = {
            'vertices': {},
            'radius': radius,
        }

        def gen_random_coordinates(v):
            r = 2*np.log(float(n) / self.graph.degree(v))
            phi = np.random.uniform(0.0, 2*np.pi)
            return (r, phi)

        def delta_phi(r, R):
            # compute \Delta\phi from equality
            # \cosh(R) = \cosh^2(r) - \sinh^2(r) * \cos(\Delta\phi)
            cos_delta_phi = (np.cosh(r)**2 - np.cosh(R)) / np.sinh(r)**2
            cos_delta_phi = max(-1., min(1., cos_delta_phi))
            return np.arccos(cos_delta_phi)

        if init_embedding is not None:
            init_embedding_fringe = fringe(self.graph, init_embedding.embedding['vertices'])
            for v in self.graph.nodes():
                if v in init_embedding.embedding['vertices']:
                    self.embedding['vertices'][v] = init_embedding.embedding['vertices'][v]
                #elif any([_n in init_embedding.embedding['vertices'] for _n in self.graph.neighbors(v)]):
                elif v in init_embedding_fringe:
                    neigh = list(set(self.graph.neighbors(v)) & set(init_embedding.embedding['vertices']))[0]
                    r_n, phi_n = init_embedding.embedding['vertices'][neigh]
                    delta_phi_n = delta_phi(r_n, radius)
                    r_random, phi_zzz = gen_random_coordinates(v)
                    phi_random_with_delta = np.random.uniform(4*np.pi + phi_n - delta_phi_n, 4*np.pi + phi_n + delta_phi_n)
                    self.embedding['vertices'][v] = (r_random, phi_random_with_delta)
                else:
                    self.embedding['vertices'][v] = gen_random_coordinates(v)
        else:
            for v in self.graph.nodes():
                self.embedding['vertices'][v] = gen_random_coordinates(v)

        self._vertex2index = {}
        for i, v in enumerate(sorted(self.embedding['vertices'])):
            self._vertex2index[v] = i

    def get_vertex_embedding(self, v):
        return self.embedding['vertices'][v]

    def __get_vertex_index(self, v):
        return self._vertex2index[v]

    def get_state_vector(self):
        x = np.zeros(self.__expected_vector_size)
        vertices = self.embedding['vertices']
        for v, i in self._vertex2index.iteritems():
            r, phi = vertices[v]
            x[2*i] = r
            x[2*i+1] = phi
        if self.fit_radius:
            x[-1] = self.embedding['radius']
        return x

    def set_state_vector(self, x, fixed_vertices=set()):
        if x.shape[0] != self.__expected_vector_size:
            raise Exception("Wrong number of elements in vector: {} instead of {}".format(x.shape[0], self.__expected_vector_size))
        vertices = self.embedding['vertices']
        for v, i in self._vertex2index.iteritems():
            if v in fixed_vertices:
                # do not update fixed vertices
                continue
            r = x[2*i]
            phi = x[2*i+1]
            vertices[v] = (r, phi)
        if self.fit_radius:
            self.embedding['radius'] = x[-1]

    def get_distance_info(self, edges_batch,fixed_vertices=set()):
        distance_info = dict()
        distance_info['radius'] = self.embedding['radius']
        distance_info['fit_radius'] = self.fit_radius
        distance_info['distances'] = dict()
        distance_info['distance_gradients'] = dict()

        for edge, w, mult in edges_batch:
            if all([v in fixed_vertices for v in edge]):
                continue
            e_v1, e_v2 = [self.get_vertex_embedding(v) for v in edge]
            distance_info['distances'][edge] = distance(e_v1, e_v2)
            # prepare sparse matrix
            dr1, dphi1, dr2, dphi2 = grad_distance(e_v1, e_v2)
            v1, v2 = edge
            x0 = np.zeros(self.__expected_vector_size)
            i_v1 = self.__get_vertex_index(v1)
            i_v2 = self.__get_vertex_index(v2)
            x0[2*i_v1] = dr1
            x0[2*i_v1+1] = dphi1
            x0[2*i_v2] = dr2
            x0[2*i_v2+1] = dphi2

            distance_info['distance_gradients'][edge] = sparse.csr_matrix(x0)

        return distance_info

