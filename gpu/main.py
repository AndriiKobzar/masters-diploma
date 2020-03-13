import numpy as np
from algorithm import *
from numba import cuda, jit

h = 0.6
alpha = 0.6
t = 1.
n = 20

points_count = 200
calculations_per_point = 100


@njit(cache=True)
def sigma(x): return np.sin(x) ** 2 + 0.05


@njit(cache=True)
def sigmaDerivative(x): return np.sin(2 * x)


@njit(cache=True)
def secondSigmaDerivative(x): return 2 * np.cos(2 * x)


@cuda.jit
def kernel():
    ty = cuda.blockIdx.x
    u = int(ty / points_count)
    return get_delta(t, n, h, alpha, sigma, sigmaDerivative, secondSigmaDerivative, u)


print(np.average([get_delta(t, n, h, alpha, sigma, sigmaDerivative, secondSigmaDerivative, 0.1) for x in range(0, 1)]))


def run_on_gpu():
    return kernel[points_count, calculations_per_point]()


#print(run_on_gpu())
