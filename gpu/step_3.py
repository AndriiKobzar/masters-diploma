import scipy.special as s
import numpy as np
import utils

def c(h):
  dividend = 2 * h * s.gamma(1.5 - h)
  divisor = s.gamma(h + 0.5) * s.gamma(2 - 2 / h)
  return (h - 0.5) * np.sqrt(dividend / divisor)

def r(h, t, N, alpha):
  delta = t / N
  def matrixElement(n, k):
    if (k>n): return 0
    return c(h)*np.power(delta, h-0.5) * utils.sumOperator(k+1, n, lambda i: np.exp(-alpha*(n-i-1)*delta)*np.power((i+1)/(k+1), h-0.5)*np.power(i-k, h-1.5))
  return np.fromfunction(np.vectorize(matrixElement), (N,N), dtype=int)
  