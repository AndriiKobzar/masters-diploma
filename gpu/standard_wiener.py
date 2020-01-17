from numpy.random import normal


def std_wiener(derivation, amount):
    return normal(0, derivation, [amount])
