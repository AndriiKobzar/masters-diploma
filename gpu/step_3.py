import scipy.special as sp
import numpy as np
from utils import gauss_kronrod_integrate
from math import exp


def c(h):
    dividend = 2 * h * sp.gamma(1.5 - h)
    divisor = sp.gamma(h + 0.5) * sp.gamma(2 - 2 / h)
    return (h - 0.5) * np.sqrt(dividend / divisor)


def kernel(h):
    return lambda s, t: (
        c(h) * (s ** (0.5 - h)) * gauss_kronrod_integrate(lambda u: (u**(h-0.5))*((u-s)**(h-1.5)), s, t) if s < t else 0
    )


def dvbyt(h, alpha):
    current_kernel = kernel(h)
    return lambda u, t: current_kernel(t, u) - alpha*exp(-alpha*t) * gauss_kronrod_integrate(lambda s: exp(alpha*s)*current_kernel(s, u), u, t)
