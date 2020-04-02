from utils import *
import step_3
from math import exp, pow


def get_delta(t, N, h, s_array, alpha, y, sbm_increments, sigma, sigmaDerivative, secondSigmaDerivative):
    eta = 1 / sum(map(lambda a: a ** 2, s_array))

    first = 2 * eta * step_3.c(h) * sum(map(lambda k: pow(k * t / N, 0.5 - h) * quad_integrate(
        lambda s: sigma(step_function(y, t / N, s)) * sigmaDerivative(step_function(y, t / N, s)) *
                  quad_integrate(
                      lambda v: integrand(t, N, h, alpha, s, v, k),
                      k * t / N, s),
        k * t / N, t) * sbm_increments[k], range(1, N - 1)))

    sndsn = gauss_kronrod_integrate(
        lambda tau: step_function(s_array, t / N, tau) * dvbeta(tau, alpha, eta, t, N, h, sigma, sigmaDerivative,
                                                                secondSigmaDerivative, s_array, y), 0, t)
    return first - sndsn


@njit(cache=True)
def integrand(t, n, h, alpha, s, v, k):
    print(k)
    return exp(-alpha * (s - v)) * pow(v, h - 0.5) * pow(v - k * t / n, h - 1.5)


def dvbeta(v, alpha, eta, t, n, h, sigma, sigma_der, sigma_second_der, s, y):
    return quad_integrate(lambda q: quad_integrate(dvbeta_integrand, 0, t, (
        q, v, eta, alpha, t, n, h, sigma, sigma_der, sigma_second_der, s, y)), 0, t)


def dvbeta_integrand(tau, q, v, eta, alpha, t, n, h, sigma, sigma_der, sigma_second_der, s, y):
    return dvbeta_integrand_optimized(tau, q, eta, t, n, sigma, sigma_der, sigma_second_der, s, y) \
           * step_3.dvbyt(h, alpha, v, tau) * step_3.dvbyt(h, alpha, q, tau)


@njit()
def dvbeta_integrand_optimized(tau, q, eta, t, n, sigma, sigma_der, sigma_second_der, s, y):
    return -4 * (eta ** 2) * step_function(s, t / n, q) * (sigma_der(step_function(y, t / n, tau)) ** 2 +
                                                           sigma(step_function(y, t / n,
                                                                               tau)) * sigma_second_der(
                step_function(y, t / n, tau)))
