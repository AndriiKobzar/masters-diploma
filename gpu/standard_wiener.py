from numba import njit
from numpy.random import normal


@njit()
def std_wiener(derivation, amount):

    return [normal(0, derivation) for x in range(0, amount)]
