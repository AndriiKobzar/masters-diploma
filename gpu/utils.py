# lower - inclusive
# higher - exclusive
def sumOperator(lower, higher, func):
  return sum(map(func, range(lower, higher)))