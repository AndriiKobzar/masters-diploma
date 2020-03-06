from math import floor

from numba import jit, cfunc
from scipy.integrate import quad

GAUSS_KRONROD_ERR_ABS = 0
GAUSS_KRONROD_ERR_REL = 0.001


def sumOperator(lower, higher, func):
    return sum(map(func, range(lower, higher)))


# def gauss_kronrod_integrate(f, lower, upper):
#     def wrapper(x, params): return f(x)
#     return qng(gsl_function(wrapper, None), lower, upper, GAUSS_KRONROD_ERR_REL, GAUSS_KRONROD_ERR_ABS)[1]
def gauss_kronrod_integrate(f, lower, upper):
    lower = lower if lower != 0 else 0.01
    # nb_integrand = cfunc("float64(float64)")(f)
    integration = quad(f, lower, upper, epsrel=GAUSS_KRONROD_ERR_REL)
    # print(integration)
    return integration[0]


def step_function(array, interval):
    def step(x):
        index = int(floor(x/interval))
        if index >= len(array):
            return array[-1]
        return array[index]
    return step
