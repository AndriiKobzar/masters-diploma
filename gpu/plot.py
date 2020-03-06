import matplotlib.pyplot as plt
import numpy as np

from step_1 import get_observations

n = 400
observations, increments = get_observations(0.6, 0.9, n, 1)


def get_fbm(increments):
    result = [0]
    value = 0
    for i in increments:
        value += i
        result.append(value)
    return result


fig, ax = plt.subplots()
x = list(np.linspace(0, 1, n))
ax.plot(x, observations, c="blue")
ax.plot(x, get_fbm(increments), c="green")
ax.set_xlabel("t")
ax.set_ylabel("Y")
plt.show()
