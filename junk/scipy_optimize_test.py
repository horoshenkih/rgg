import numpy as np
from scipy.optimize import minimize

class F:
    def __init__(self, bias):
        self.bias = bias
    def __call__(self, x):
        return sum([(xi - bi)**2 / 2. for (bi, xi) in zip(self.bias, x)])

class GradF(F):
    def __call__(self, x):
        return np.array([(xi - bi) for (bi, xi) in zip(self.bias, x)])

bias = (1., 2., 3.)
f = F(bias)
grad_f = GradF(bias)

x0 = (1., 1., 1.)

print minimize(f, x0, jac=grad_f)
