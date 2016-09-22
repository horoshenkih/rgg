from abc import ABCMeta, abstractmethod

import numpy as np

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
    def __init__(self, graph, radius=None):
        self.graph = graph
        n = self.graph.number_of_nodes()
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

    def get_state_vector(self):
        x = []
        for v in sorted(self.embedding['vertices'].keys()):
            r, phi = self.embedding['vertices'][v]
            x.append(r)
            x.append(phi)

        return np.array(x)

    def set_state_vector(self, x):
        if x.shape[0] != self.__expected_vector_size:
            raise Exception("Wrong number of elements in vector: {} instead of {}".format(x.shape[0], self.__expected_vector_size))
        for i, v in enumerate(self.embedding['vertices'].keys()):
            r = x[2*i]
            phi = x[2*i+1]
            self.embedding['vertices'][v] = (r, phi)

    def get_distance_info(self, edges_batch):
        distance_info = dict()
        distance_info['radius'] = self.embedding['radius']
        distance_info['distances'] = dict()
        distance_info['distance_gradients'] = dict()

        for edge in edges_batch():
            distance_info['distances'][edge] = distance(*edge)
            distance_info['distance_gradients'][edge] = grad_distance(*edge)

        return distance_info

