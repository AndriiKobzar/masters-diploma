import utils
import step_3
import numpy as np

def get_delta(t, N, h, s, alpha, y, r, fbmIncrements, sigma, sigmaDerivative, secondSigmaDerivative):
  delta = t/N
  eta=1/(delta*utils.sumOperator(0, N, lambda i: s[i]**2))
  negativeSumSnDsn = 16 * eta**2 * delta**3 * utils.sumOperator(0, N, lambda k1: \
      utils.sumOperator(0,N, lambda k2: \
        utils.sumOperator(0, N, lambda k3: \
          utils.sumOperator(0, N, lambda i: sigma(y[i])*sigmaDerivative(y[i])*r[k2][i]) * \
          utils.sumOperator(0, N, lambda j: sigma(y[j])*sigmaDerivative(y[j])*r[k1][j]) * \
          (sigmaDerivative(y[k3])**2 + secondSigmaDerivative(y[k3])*sigma(y[k3])) * \
          r[k1][k3]*r[k2][k3] \
          )))
  return negativeSumSnDsn + utils.sumOperator(1, N, lambda n: \
    utils.sumOperator(0, n-1, lambda k: \
      utils.sumOperator(k+1, n+1, lambda i: 2*eta*step_3.c(h)*delta**(h+0.5)*fbmIncrements[k]*np.exp(-alpha*delta*(n-i-1))*np.power((i+1)/(k+1),h-0.5)*np.power(i-k, h-1.5))))