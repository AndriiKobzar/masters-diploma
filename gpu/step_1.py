import fbm_increments as fbmi

def get_observations(alpha, h, n, t):
  fbm_increments = fbmi.get_fbm_increments(t, n, h)
  increment = t / n
  currentValue = 0
  currentCoordinate = 0
  result = [0]*n

  for index in range(1, n):
    currentValue = currentValue - alpha * currentValue * increment + fbm_increments[index-1]

    currentCoordinate += increment
    result[index] = currentCoordinate
  return (result, fbm_increments)
