import utils
from step_3 import dvbyt


def s_array(t, n, sigma, sigmaDerivative, y, alpha, h):
    integrand = dvb_sigma_square_y(t, n, sigma, sigmaDerivative, y, alpha, h)

    def array_elem(k):
        return 2 * utils.gauss_kronrod_integrate(lambda x: integrand(k, x), 0, t)

    return list(map(array_elem, range(0, n)))


def dvb_sigma_square_y(t, n, sigma, sigma_der, y, alpha, h):
    step_y = utils.step_function(y, t / n)
    current_dvbyt = dvbyt(h, alpha)
    return lambda k, s: sigma(step_y(s)) * sigma_der(step_y(s)) * current_dvbyt(k * t / n, s)
