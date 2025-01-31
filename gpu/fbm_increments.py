import numpy as np
from numba import jit, njit


def get_lambda(h, n):
    m = 2 * n - 2
    c = [0] * m
    g = 2 * h

    for i in range(0, n):
        c[i] = fbc(g, i)
    for i in range(1, n):
        c[-i] = c[i]
    return np.sqrt(np.fft.fft(c).real)


def fbc(g, x):
    return ((x + 1) ** g + abs(x - 1) ** g - 2 * x ** g) / 2


def get_fgn(lambda_array):
    m = len(lambda_array)
    inverse_transform_on_random = np.fft.ifft(np.random.normal(size=m))
    a = np.array(inverse_transform_on_random) * np.array(lambda_array)
    result = np.real(np.fft.fft(a))
    return result[:m // 2]


def get_fbm_increments(t, n, h):
    lambda_array = get_lambda(h, n)
    fgn = get_fgn(lambda_array)
    return fgn * ((t / n) ** h)
