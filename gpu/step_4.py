import utils
import numpy as np

def s(t, n, sigma, sigmaDerivative, observations, r):
  delta = t/n

  def array_elem(k):
    return 2 * delta * utils.sumOperator(0, n,
     lambda i: sigma(i)*sigmaDerivative(i)*r.item((i,int(k))))
  return np.fromfunction(np.vectorize(array_elem), [n])
