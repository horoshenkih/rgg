from abc import ABCMeta, abstractmethod

import numpy as np

class EmbeddingModel:
    __metaclass__ = ABCMeta
    @abstractmethod
    def __init__(self):
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
    def get_distance_info(self, edges):
        """Get distances, gradients and other information for given set of edges"""
        pass

class PoincareFixedRadius(EmbeddingModel):
    def __init__(self):
        print 'OK'
