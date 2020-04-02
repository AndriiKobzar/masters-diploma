from numba import njit, jit

import utils
from step_3 import dvbyt


def s_array(t, n, sigma, sigma_der, y, alpha, h):

    def array_elem(k):
        print("s array ", k)
        return 2 * utils.quad_integrate(dvb_sigma_square_y, 0, t, (t, n, sigma, sigma_der, y, alpha, h, k))

    return list(map(array_elem, range(0, n)))


def dvb_sigma_square_y(s, t, n, sigma, sigma_der, y, alpha, h, k):
    return sigma(utils.step_function(y, t / n, s)) * sigma_der(utils.step_function(y, t / n, s)) * dvbyt(h, alpha, k * t / n, s)
