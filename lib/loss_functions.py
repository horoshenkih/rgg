from abc import ABCMeta, abstractmethod

import numpy as np

class LossFunction:
    __metaclass__ = ABCMeta
    @abstractmethod
    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def loss(self, edges_batch, distance_info):
        pass

    @abstractmethod
    def loss_gradient(self, edges_batch, distance_info):
        pass

class SmoothEdgePredictor:
    def predict(self, d, r, beta):
        return 1. / (1 + np.exp((d - r)*float(beta)))

    def prediction_gradient_d(self, d, r, beta):
        p = self.predict(d, r, beta)
        return -np.exp((d - r)*float(beta)) * beta * p**2

    def prediction_gradient_r(self, d, r, beta):
        return -self.prediction_gradient_d(d, r, beta)

class MSE(LossFunction):
    def __init__(self, binary_edges=False, beta=1.):
        self.binary_edges = binary_edges
        self.edge_predictor = SmoothEdgePredictor()
        self.beta = beta

    def loss(self, edges_batch, distance_info):
        r = distance_info['radius']
        total_loss = 0.
        for e, w, mult in edges_batch:
            if self.binary_edges:
                w = 1. if w else 0.
            d = distance_info['distances'][e]
            p = self.edge_predictor.predict(d, r, self.beta)
            total_loss += (p - w)**2 * mult
        return total_loss

    def loss_gradient(self, edges_batch, distance_info):
        r = distance_info['radius']
        grad = None
        for e, w, mult in edges_batch:
            if self.binary_edges:
                w = 1. if w else 0.
            d = distance_info['distances'][e]
            p = self.edge_predictor.predict(d, r, self.beta)
            dp = self.edge_predictor.prediction_gradient_d(d, r, self.beta)
            dd = distance_info['distance_gradients'][e].toarray()[0]
            d_grad = 2. * dd * dp * (p - w) * mult
            if distance_info['fit_radius']:
                d_grad[-1] = 2 * self.edge_predictor.prediction_gradient_r(d, r, self.beta) * (p - w) * mult

            if grad is None:
                grad = d_grad
            else:
                grad += d_grad
        return grad

