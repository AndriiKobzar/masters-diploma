from pygsl.integrate import qng
from math import floor

GAUSS_KRONROD_ERR_ABS = 0.001
GAUSS_KRONROD_ERR_REL = 0.001


def sumOperator(lower, higher, func):
    return sum(map(func, range(lower, higher)))


def gauss_kronrod_integrate(f, lower, upper):
    return qng(f, lower, upper, GAUSS_KRONROD_ERR_REL, GAUSS_KRONROD_ERR_ABS)


def step_function(array, interval):
    return lambda x: array[int(floor(x/interval))]
