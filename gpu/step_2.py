def integral_of_square_sigma(t, sigma, observations):
  delta = t/len(observations)
  return sum(map(lambda x: delta*sigma(x)**2, observations))