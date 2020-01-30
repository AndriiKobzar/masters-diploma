import fbm_increments as fbmi


def get_observations(alpha, h, n, t):
    fbm_increments = fbmi.get_fbm_increments(t, n, h)
    increment = t / n
    currentValue = 0
    result = [0] * n

    for index in range(1, n):
        currentValue = currentValue - alpha * currentValue * increment + fbm_increments[index - 1]
        result[index] = currentValue
    return (result, fbm_increments)
