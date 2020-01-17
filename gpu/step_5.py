import utils
import step_3
from math import exp, pow


def get_delta(t, N, h, s_array, alpha, y, sbm_increments, sigma, sigmaDerivative, secondSigmaDerivative):
    eta = 1 / sum(map(lambda a: a ** 2, s_array))
    inner_integrand = integrand(t, N, h, alpha)
    step_y = utils.step_function(y, t / N)
    first = 2 * eta * step_3.c(h) * sum(map(lambda k: pow(k * t / N, 0.5 - h) *
                                                      utils.gauss_kronrod_integrate(
                                                          lambda s: sigma(step_y(s)) * sigmaDerivative(step_y(s)) *
                                                                    utils.gauss_kronrod_integrate(lambda v: inner_integrand(s, v, k),
                                                                              k * t / N, s),
                                                          k * t / N, t) * sbm_increments[k], range(1, N - 1)))
    step_s = utils.step_function(s_array, t / N)

    def dvbeta(v):
        dvbyt = step_3.dvbyt(h, alpha)
        return -4 * eta**2 * utils.gauss_kronrod_integrate(lambda q: step_s(q)*utils.gauss_kronrod_integrate(lambda tau: (sigmaDerivative(step_y(tau))**2 +
                                                                                sigma(step_y(tau))*secondSigmaDerivative(step_y(tau)))
                                                                     * dvbyt(v, tau)*dvbyt(q,tau), 0, t), 0, t)

    sndsn = utils.gauss_kronrod_integrate(lambda tau: step_s(tau)*dvbeta(tau), 0, t)
    return first - sndsn


def integrand(t, n, h, alpha):
    return lambda s, v, k: exp(-alpha * (s - v)) * pow(v, h - 0.5) * pow(v - k * t / n, h - 1.5)
