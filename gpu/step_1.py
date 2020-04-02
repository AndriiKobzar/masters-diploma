from numba import njit, jit
from fbm_increments import get_fbm_increments



def get_observations(alpha, h, n, t):
    fbm_increments = get_fbm_increments(t, n, h)
    increment = t / n
    currentValue = .0
    result = [0.] * n

    for index in range(1, n):
        currentValue = currentValue - alpha * currentValue * increment + fbm_increments[index - 1]
        result[index] = currentValue
    return result
