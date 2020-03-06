import step_1
import step_2
import step_4
import step_5
from standard_wiener import *


def get_delta(t, n, h, alpha, sigma, sigmaDerivative, secondSigmaDerivative, u):
    observations = step_1.get_observations(alpha, h, n, t)
    print("observations")
    integral = step_2.integral_of_square_sigma(t, sigma, observations)
    if integral <= u:
        return 0
    s = step_4.s_array(t, n, sigma, sigmaDerivative, observations, alpha, h)
    print("array s")
    sbm_increments = std_wiener(t / n, n)
    delta = step_5.get_delta(t, n, h, s, alpha, observations, sbm_increments, sigma, sigmaDerivative,
                             secondSigmaDerivative)
    print(delta)
    return delta
