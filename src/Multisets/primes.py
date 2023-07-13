import itertools

def eratosthenes_enum():
  d = dict()
  yield 2
  for q in itertools.islice(itertools.count(3), 0, None, 2):
    p = d.pop(q, None)
    if p is None:
      d[q*q] = q
      yield q
    else:
      x = p + q
      while x in d or not (x&1):
        x += p
      d[x] = p

def next_larger_prime(lowerbound):
  for prime in eratosthenes_enum():
    if prime > lowerbound:
      return prime
