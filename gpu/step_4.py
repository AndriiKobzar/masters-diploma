import utils
from step_3 import dvbyt


def s_array(t, n, sigma, sigma_der, y, alpha, h):

    def array_elem(k):
        print("s array ", k)

        def integrand(x):
            return dvb_sigma_square_y(t, n, sigma, sigma_der, y, alpha, h, k, x)

        return 2 * utils.gauss_kronrod_integrate(integrand, 0, t)

    return list(map(array_elem, range(0, n)))


def dvb_sigma_square_y(t, n, sigma, sigma_der, y, alpha, h, k, s):
    step_y = utils.step_function(y, t / n)
    return sigma(step_y(s)) * sigma_der(step_y(s)) * dvbyt(h, alpha, k * t / n, s)
