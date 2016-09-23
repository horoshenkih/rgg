from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import sparse

from graph import distance, grad_distance

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
    def __init__(self, graph, radius=None, fit_radius=False):
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
        for v in self.graph.nodes():
            r = 2*np.log(float(n) / self.graph.degree(v))
            phi = np.random.uniform(0.0, 2*np.pi)
            self.embedding['vertices'][v] = (r, phi)
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

    def set_state_vector(self, x):
        if x.shape[0] != self.__expected_vector_size:
            raise Exception("Wrong number of elements in vector: {} instead of {}".format(x.shape[0], self.__expected_vector_size))
        vertices = self.embedding['vertices']
        for v, i in self._vertex2index.iteritems():
            r = x[2*i]
            phi = x[2*i+1]
            vertices[v] = (r, phi)
        if self.fit_radius:
            self.embedding['radius'] = x[-1]

    def get_distance_info(self, edges_batch):
        distance_info = dict()
        distance_info['radius'] = self.embedding['radius']
        distance_info['fit_radius'] = self.fit_radius
        distance_info['distances'] = dict()
        distance_info['distance_gradients'] = dict()

        for edge, w, mult in edges_batch:
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

