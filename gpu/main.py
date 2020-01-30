import numpy as np
import step_1
import step_2
import step_4
import step_5
from standard_wiener import *
import matplotlib.pyplot as plt

h = 0.6
alpha = 0.6
t = 1
n = 2000


def sigma(x): return np.sin(x) ** 2 + 0.05


def sigmaDerivative(x): return np.sin(2 * x)


def secondSigmaDerivative(x): return 2 * np.cos(2 * x)


def get_delta(u):
    (observations, fbmIncrements) = step_1.get_observations(alpha, h, n, t)
    integral = step_2.integral_of_square_sigma(t, sigma, observations)
    if integral <= u:
        return 0
    s = step_4.s_array(t, n, sigma, sigmaDerivative, observations, alpha, h)
    sbm_increments = std_wiener(t / n, n)
    delta = step_5.get_delta(t, n, h, s, alpha, observations, sbm_increments, sigma, sigmaDerivative,
                             secondSigmaDerivative)
    print(delta)
    return delta


def plot():
    (observations, fbmIncrements) = step_1.get_observations(alpha, h, n, t)
    fig, ax = plt.subplots()
    ax.plot(observations)
    plt.show()


print(np.average([get_delta(0.1) for x in range(0, 100)]))

#plot()
