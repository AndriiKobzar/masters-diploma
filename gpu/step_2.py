def integral_of_square_sigma(t, sigma, observations):
    result = 0
    for y in observations:
        result += sigma(y) ** 2
    return result
