from pygsl.integrate import qng
from pygsl.gsl_function import gsl_function
from math import floor

GAUSS_KRONROD_ERR_ABS = 0
GAUSS_KRONROD_ERR_REL = 0.001


def sumOperator(lower, higher, func):
    return sum(map(func, range(lower, higher)))


def gauss_kronrod_integrate(f, lower, upper):
    def wrapper(x, params): return f(x)
    return qng(gsl_function(wrapper, None), lower, upper, GAUSS_KRONROD_ERR_REL, GAUSS_KRONROD_ERR_ABS)[1]


def step_function(array, interval):
    return lambda x: array[int(floor(x/interval))]
