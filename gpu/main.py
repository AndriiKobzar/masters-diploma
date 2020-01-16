import numpy as np
import step_1
import step_2
import step_3
import step_4
import step_5

h=0.6
alpha=0.5
t=1
n=10
u=1

def sigma(x): return np.sin(x)**2 + 0.05
def sigmaDerivative(x): return np.sin(2*x)
def secondSigmaDerivative(x): return 2 * np.cos(2*x)

def getDelta():
  (observations, fbmIncrements) = step_1.get_observations(alpha, h, n, t)
  integral = step_2.integral_of_square_sigma(t, sigma, observations)
  while integral>u:
    (observations, fbmIncrements) = step_1.get_observations(alpha, h, n, t)
    integral = step_2.integral_of_square_sigma(t, sigma, observations)
  print(observations)
  r = step_3.r(h, t, n, alpha)
  print(r)
  s = step_4.s(t,n, sigma, sigmaDerivative, observations, r)
  print(s)
  return step_5.get_delta(t, n, h, s, alpha, observations, r, fbmIncrements, sigma, sigmaDerivative, secondSigmaDerivative)
  

print(getDelta())