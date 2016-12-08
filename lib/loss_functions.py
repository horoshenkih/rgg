from abc import ABCMeta, abstractmethod

import numpy as np


class LossFunction:
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, **kwargs):
        pass

    @abstractmethod
    def loss(self, edges_batch, distance_info):
        return 0.

    @abstractmethod
    def loss_gradient(self, edges_batch, distance_info):
        pass


class SmoothEdgePredictor:
    def predict(self, d, r, beta):
        return 1. / (1 + np.exp((d - r) * float(beta)))

    def prediction_gradient_d(self, d, r, beta):
        p = self.predict(d, r, beta)
        return -np.exp((d - r) * float(beta)) * beta * p ** 2

    def prediction_gradient_r(self, d, r, beta):
        return -self.prediction_gradient_d(d, r, beta)


class SmoothEdgeLoss(LossFunction):
    __metaclass__ = ABCMeta

    def __init__(self, binary_edges=False, beta=1.):
        self.binary_edges = binary_edges
        self.edge_predictor = SmoothEdgePredictor()
        self.beta = beta
        self.EPS = 1e-15

    @abstractmethod
    def elementary_loss(self, true, predicted):
        """
        Loss in individual object (e.g. (true - predicted)**2 for MSE)
        :param true: true value
        :param predicted: predicted value
        :return: value of loss
        """
        return 0.

    @abstractmethod
    def elementary_loss_gradient(self, true, predicted):
        """
        Gradient of loss wrt predicted value (e.g. 2. * (predicted - true) for MSE)
        :param true: true value
        :param predicted: predicted value
        :return: value of gradient
        """
        return 0.

    def loss(self, edges_batch, distance_info):
        r = distance_info['radius']
        total_loss = 0.
        for e, w, mult in edges_batch:
            if self.binary_edges:
                w = 1. if w else 0.
            d = distance_info['distances'][e]
            p = self.edge_predictor.predict(d, r, self.beta)
            total_loss += self.elementary_loss(w, p) * mult
        return total_loss

    def loss_gradient(self, edges_batch, distance_info):
        r = distance_info['radius']
        grad = None
        for e, w, mult in edges_batch:
            if self.binary_edges:
                w = 1. if w else 0.
            if e not in distance_info['distances']:
                continue
            d = distance_info['distances'][e]
            p = self.edge_predictor.predict(d, r, self.beta)
            dp = self.edge_predictor.prediction_gradient_d(d, r, self.beta)
            dd = distance_info['distance_gradients'][e].toarray()[0]
            d_grad = dd * dp * mult * self.elementary_loss_gradient(w, p)
            if distance_info['fit_radius']:
                d_grad[-1] = self.edge_predictor.prediction_gradient_r(d, r, self.beta) * mult * self.elementary_loss_gradient(w, p)

            if grad is None:
                grad = d_grad
            else:
                grad += d_grad
        return grad


class MSE(SmoothEdgeLoss):
    def elementary_loss(self, true, predicted):
        return (predicted - true)**2

    def elementary_loss_gradient(self, true, predicted):
        return 2.*(predicted - true)


class LogLoss(SmoothEdgeLoss):
    def elementary_loss(self, true, predicted):
        p = max(self.EPS, predicted)
        p = min(p, 1-self.EPS)
        return - true * np.log(p) - (1. - true) * np.log(1 - p)

    def elementary_loss_gradient(self, true, predicted):
        p = max(self.EPS, predicted)
        p = min(p, 1 - self.EPS)
        return (p - true) / p / (1-p)
