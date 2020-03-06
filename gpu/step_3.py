import numpy as np
from utils import gauss_kronrod_integrate
from math import exp, gamma


def c(h):
    dividend = 2 * h * gamma(1.5 - h)
    divisor = gamma(h + 0.5) * gamma(2 - 2 / h)
    return (h - 0.5) * np.sqrt(dividend / divisor)


def kernel(h):
    return lambda s, t: (
        c(h) * (s ** (0.5 - h)) * gauss_kronrod_integrate(
            lambda u: (u ** (h - 0.5)) * ((u - s if u - s == 0 else 0.01) ** (h - 1.5)), s, t) if s < t else 0
    )


def dvbyt(h, alpha, u, t):
    return gauss_kronrod_integrate(lambda s: dvbyt_integrand(h, alpha, u, t, s), u, t) if u < t else 0


def dvbyt_integrand(h, alpha, u, t, s):
    return c(h) * exp(-alpha * t) * (close_to_zero(u) ** (0.5 - h)) * ((close_to_zero(u - s)) ** (h - 0.5)) * exp(alpha * s) * (close_to_zero(s - u) ** (h - 1.5))


def close_to_zero(x):
    return 0.01 if x == 0 else x
