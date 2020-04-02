def dvbeta(v):
    return -4 * (eta ** 2) * utils.quad_integrate(
        lambda q:  utils.quad_integrate(lambda tau: step_s(q) * (sigmaDerivative(step_y(tau)) ** 2 +
                                                                 sigma(step_y(tau)) * secondSigmaDerivative(
                    step_y(tau))) * step_3.dvbyt(h, alpha, v, tau) * step_3.dvbyt(h, alpha, q, tau), 0, t), 0, t)


sndsn = utils.quad_integrate(lambda tau: step_s(tau) * dvbeta(tau), 0, t)