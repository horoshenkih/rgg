from abc import ABCMeta, abstractmethod

import numpy as np

class LossFunction:
    __metaclass__ = ABCMeta
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def generate_edges_set(self):
        pass

    @abstractmethod
    def update_state_vector(self, edges, x, distance_info):
        """Return new state vector"""
        pass
