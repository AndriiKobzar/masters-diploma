import numpy as np

def get_lambda(h, n):
  m = 2*n-2
  c = [0]*m
  g = 2*h

  def fbc(x): return ((x+1)**g + abs(x-1)**g - 2*x**g)/2
  for x in range(0, n):
    c[x] = fbc(x)
  for x in range(1, n):
    c[-x] = c[x]
  return np.sqrt(np.fft.fft(c).real)

def get_fgn(lambda_array):
  m=len(lambda_array)
  inverse_transform_on_random = np.fft.ifft(np.random.normal(size=m))
  a = np.array(inverse_transform_on_random)*np.array(lambda_array)
  result = np.real(np.fft.fft(a))
  return result[:m//2]

def get_fbm_increments(t, n, h):
  lambda_array = get_lambda(h, n)
  fgn = get_fgn(lambda_array)
  return fgn*((t/n)**h)
