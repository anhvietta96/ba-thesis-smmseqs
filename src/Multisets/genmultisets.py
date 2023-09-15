#!/usr/bin/env python3

import sys, argparse, textwrap
from itertools import permutations
from math import factorial
from primes import next_larger_prime

def binom_evaluate(n,k):
  if k == 0 or n == k:
    return 1
  if n == 0:
    return 0
  result = 1
  for j in range(1,k+1):
    result = (result * (n - k + j))//j
  return result

# determine number of multisets with m elements chosen from a set of size k
# the formula which is from a wikipedia entry about multisets

def multisets_number(m,k):
  return binom_evaluate(k + m - 1, m)

# function to recursively compute the number of multisets of size m over
# the alphabet 0, 1, ..., k-1, where k is the number of symbols

def multisets_count(multiset_size,num_symbols):
  def count_rec(m,k):
    if m == 0:
      return 1
    if k == 0:
      return 0
    count = 0
    for i in range(0,m+1):
      count += count_rec(m-i,k-1)
    return count
  return count_rec(multiset_size,num_symbols)

# generate multisets and yield them one after the other

def multisets_enumerate_stack(multiset_size,num_symbols):
  multiset = [0] * multiset_size
  stack = [(multiset_size,num_symbols)]
  while stack:
    (m,k) = stack.pop()
    if m == 0:
      yield multiset
    elif k > 0:
      stack.append((m,k-1))
      for i in range(1,m+1):
        stack.append((m-i,k-1))
        for d in range(i):
          assert d < m
          multiset[multiset_size - m + d] = num_symbols - k

def multisets_enumerate_loops(multiset_size,num_symbols):
  if multiset_size == 2:
    for a in range(0,num_symbols):
      for b in range(a,num_symbols):
        yield a, b
  elif multiset_size == 3:
    for a in range(0,num_symbols):
      for b in range(a,num_symbols):
        for c in range(b,num_symbols):
          yield a, b, c
  elif multiset_size == 4:
    for a in range(0,num_symbols):
      for b in range(a,num_symbols):
        for c in range(b,num_symbols):
          for d in range(c,num_symbols):
            yield a, b, c, d
  elif multiset_size == 5:
    for a in range(0,num_symbols):
      for b in range(a,num_symbols):
        for c in range(b,num_symbols):
          for d in range(c,num_symbols):
            for e in range(d,num_symbols):
              yield a, b, c, d, e
  elif multiset_size == 6:
    for a in range(0,num_symbols):
      for b in range(a,num_symbols):
        for c in range(b,num_symbols):
          for d in range(c,num_symbols):
            for e in range(d,num_symbols):
              for f in range(e,num_symbols):
                yield a, b, c, d, e, f
  elif multiset_size == 7:
    for a in range(0,num_symbols):
      for b in range(a,num_symbols):
        for c in range(b,num_symbols):
          for d in range(c,num_symbols):
            for e in range(d,num_symbols):
              for f in range(e,num_symbols):
                for g in range(f,num_symbols):
                  yield a, b, c, d, e, f, g

def multisets_enumerate(multiset_size,num_symbols):
  if multiset_size <= 7:
    for ms in multisets_enumerate_loops(multiset_size,num_symbols):
      yield ms
  else:
    for ms in multisets_enumerate_stack(multiset_size,num_symbols):
      yield ms

def multiset_symbols(multiset):
  previous = None
  symbols = list()
  for sym in multiset:
    if previous is None or previous != sym:
      symbols.append(sym)
    previous = sym
  return symbols

def multiset_enumerate_counts(multiset):
  previous = None
  count = 0
  for sym in multiset:
    if previous is None or previous == sym:
      count += 1
    else:
      yield count
      count = 1
    previous = sym
  yield count

def multiset_is_it(intlist):
  for idx in range(1,len(intlist)):
    if intlist[idx-1] > intlist[idx]:
      return False
  return True

def partialsum(l):
  psums = [None] * len(l)
  psums[0] = l[0]
  for i in range(1,len(l)):
    psums[i] = psums[i-1] + l[i]
  return psums

def first_factors_new(k):
  previous_sum = 0
  sums = list()
  for i in range(1,k):
    previous_sum += i
    sums.append(previous_sum)
  sums.append(0)
  sums = list(reversed(sums))
  psums = partialsum(sums)
  return psums

def second_factors_new(k):
  sums = list()
  for i in range(0,k):
    sums.append((i * (2 * k - i - 1))//2)
  return sums

def multiset_hash_encode(ms,alphasize,first_factors,second_factors):
  if len(ms) == 2:
    # first = ms[0] * alphasize - (ms[0] * (ms[0]+1))//2
    return first_factors[ms[0]] + ms[1]
  elif len(ms) == 3:
    # second = ms[1] * alphasize - (ms[1] * (ms[1]+1))//2
    return first_factors[ms[0]] + \
           second_factors[ms[1]] + \
           ms[2]

def diff_func(k,r):
  l = list()
  l.append(0)
  for i in range(1,r+1):
    l.append(((k-i) * (k-i+1))//2)
  return sum(l)//2

def multiset_weights(alphasize,multiset_size):
  weights = list(range(0,alphasize))
  weights_table = [weights]
  for m in range(2,multiset_size+1):
    idx = 0
    weights = [None] * alphasize
    for multiset in multisets_enumerate(m,alphasize):
      suffixcode = 0
      for pos in range(0,m-1):
        suffixcode += weights_table[pos][multiset[m-1-pos]]
      print("{} {} {}".format(multiset,suffixcode,idx))
      assert suffixcode <= idx
      difference = idx - suffixcode
      if weights[multiset[0]] is None:
        weights[multiset[0]] = difference
      else:
        assert weights[multiset[0]] == difference
      idx += 1
    weights_table.append(weights)
  return weights_table

def multiset_hash_encode_all(m,k):
  second_factors = second_factors_new(k)
  idx = 0
  if m == 2:
    for multiset in multisets_enumerate(m,k):
      code = multiset_hash_encode(multiset,k,second_factors,None)
      assert idx == code
      print("{}\t{}".format(idx,multiset))
      idx += 1
  elif m == 3:
    first_factors = first_factors_new(k)
    print("j\tfirst\tsecond");
    for j in range(0,k):
      print("{}\t{}\t{}".format(j,first_factors[j],second_factors[j]))
    for multiset in multisets_enumerate(m,k):
      code = multiset_hash_encode(multiset,k,first_factors,second_factors)
      assert code == idx
      print("{}\t{}".format(idx,multiset))
      idx += 1
  elif m == 4:
    first_factors = first_factors_new(k)
    diff4 = dict()
    for multiset in multisets_enumerate(m,k):
      suffixcode = first_factors[multiset[1]] + \
                   second_factors[multiset[2]] + \
                   multiset[3]
      assert suffixcode <= idx
      difference = idx - suffixcode
      if multiset[0] in diff4:
        assert diff4[multiset[0]] == difference
      else:
        diff4[multiset[0]] = difference
      idx += 1
    for idx in range(0,k):
      print("factor4[{}]\t{}".format(idx,diff4[idx]))

def multiset_number_of_permutations(multiset):
  n = len(multiset)
  sum_counts = 0
  factorial_product = 1
  count_list = list()
  for count in multiset_enumerate_counts(multiset):
    factorial_product *= factorial(count)
    sum_counts += count
    count_list.append(count)
  assert sum_counts == n
  return factorial(n)//factorial_product

# enumerate all unique permutations of a multiset, using a brute force method
# this first enumerate all permutations including (possible) duplicates.
# The permutations are stored in list, which is sorted, before
# all duplicates are removed.

def multiset_permutations_brute_force(multiset):
  all_perms = list()
  for perm in permutations(multiset):
    all_perms.append(perm)
  all_perms = sorted(all_perms)
  previous_perm = None
  for perm in all_perms:
    if previous_perm is None or previous_perm != perm:
      yield list(perm)
    previous_perm = perm

# directly enumerate all unique permutations of a multiset using the
# Cool-lex algorithm of Williams, SODA 2009.

def multiset_cool_lex_enumerate_permutations(multiset):
  def linked_list2array(value_tab,head,next_tab):
    perm = list()
    while head is not None:
      perm.append(value_tab[head])
      head = next_tab[head]
    return perm
  if len(multiset_symbols(multiset)) == 1:
    yield multiset
    return
  m = len(multiset)
  assert m >= 2
  next_tab = [None] * m
  value_tab = [None] * m
  for idx in range(m):
    assert idx == 0 or multiset[idx-1] <= multiset[idx]
    value_tab[m-1-idx] = multiset[idx]
  for idx in range(m-1):
    next_tab[idx] = idx + 1
  head = 0
  yield linked_list2array(value_tab,head,next_tab)
  i = m - 2
  afteri = m - 1
  while (next_tab[afteri] is not None) or value_tab[afteri] < value_tab[head]:
    if (next_tab[afteri] is not None) and \
       value_tab[i] >= value_tab[next_tab[afteri]]:
      before_k = afteri
    else:
      before_k = i
    k = next_tab[before_k]
    next_tab[before_k] = next_tab[k]
    next_tab[k] = head
    if value_tab[k] < value_tab[head]:
      i = k
    assert i is not None
    afteri = next_tab[i]
    assert afteri is not None
    head = k
    assert head is not None
    yield linked_list2array(value_tab,head,next_tab)

def multisets_cmp(ms1,ms2):
  assert len(ms1) == len(ms2)
  for a, b in zip(ms1,ms2):
    if a < b:
      return -1
    if a > b:
      return 1
  return 0

def multiset_output(m,k,showms,showmsperm):
  count_sequences = k ** m
  count1 = multisets_number(m,k)
  count_unique = 0
  if showms or showmsperm:
    count2 = multisets_count(m,k)
    assert count1 == count2
    previous_multiset = None
    for multiset in multisets_enumerate(m,k):
      symbols = multiset_symbols(multiset)
      print("{}\tdifferent symbols: {}".format(multiset,len(symbols)))
      if previous_multiset is not None and\
         multiset_cmp(previous_multiset,multiset) >= 0:
        sys.stderr.write("unexpected order{} >= {}"
                          .format(previous_multiset,multiset))
        exit(1)
      if showmsperm:
        local_count = 0
        for perm in multiset_cool_lex_enumerate_permutations(multiset):
          local_count += 1
          print('\t{}'.format(perm))
        print("\tpermutations: {}".format(local_count))
        assert local_count == multiset_number_of_permutations(multiset)
        count_unique += local_count
    assert (not showmsperm) or count_unique == count_sequences
  print('{}\t{}\t{}\t{}\t{:.2f}%'.format(m,k,count1,count_sequences,
                                         100.0 * count1/count_sequences))

def multiset_permutations_verify(multiset):
  print("check permutation of multiset {}".format(multiset))
  perms1 = list()
  for perm in multiset_permutations_brute_force(multiset):
    perms1.append(perm)
  perms2 = list()
  for perm in multiset_cool_lex_enumerate_permutations(multiset):
    perms2.append(perm)
  perms1 = sorted(perms1)
  perms2 = sorted(perms2)
  if len(perms1) != len(perms2):
    sys.stderr.write('{}: for {}: perms of different length {} != {}\n'
                      .format(sys.argv[0],multiset,len(perms1),len(perms2)))
    sys.stderr.write('perms1={}\n'.format(perms1))
    sys.stderr.write('perms2={}\n'.format(perms2))
    exit(1)
  for idx in range(len(perms1)):
    if perms1[idx] != perms2[idx]:
      sys.stderr.write('{}: for {}: perms[{}] are different {} != {}\n'
                        .format(sys.argv[0],multiset,idx,
                                perms1[idx],perms2[idx]))
      sys.stderr.write('perms1={}\n'.format(perms1))
      sys.stderr.write('perms2={}\n'.format(perms2))
      exit(1)

def multisets_verify(alphasize):
  if alphasize is None:
    for k in range(2,20+1):
      multiset_hash_encode_all(2,k)
      multiset_hash_encode_all(3,k)
  else:
    multiset_hash_encode_all(2,alphasize)
    multiset_hash_encode_all(3,alphasize)
  exit(0)
  check_generate = 0
  for k in range(1,10+1):
    for m in range(1,8+1):
      value1 = multisets_number(m,k)
      value2 = multisets_count(m,k)
      value3 = 0
      for _ in multisets_enumerate(m,k):
        value3 += 1
      assert value1 == value2
      assert value2 == value3
      print('{}\t{}\t{}'.format(m,k,value1))
      check_generate += 1
  check_permutation = 0
  for k in range(2,6+1):
    for m in range(2,7+1):
      print("{} {}".format(m,k))
      for multiset in multisets_enumerate(m,k):
        check_permutation += 1
        multiset_permutations_verify(multiset)
  print("successfully checked generation of {} multisets"
         .format(check_generate))
  print("successfully checked {} multiset permutations"
         .format(check_permutation))

def parse_command_line():
  p = argparse.ArgumentParser(
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description="verify generators for multisets of m elements over alphabets of size k, verify the methods for enumerating permutations of multisets; output stattistics of sizes of depending on m k")
  processgroup = p.add_mutually_exclusive_group(required=False)
  processgroup.add_argument("-s","--show",help="show multisets",
                 action='store_true',default=False)
  processgroup.add_argument("-v","--verify",
                 help="verify that recursively determined numbers are correct",
                 action='store_true',default=False)
  processgroup.add_argument("-p","--perm",
                            help="show multisets and their permutations",
                            action='store_true',default=False)
  processgroup.add_argument("-w","--weights_table",
                 help="output weights table for give",
                 action='store_true',default=False)
  p.add_argument("-k","--num_symbols",
                 help="specify size of base set ",type=int,default=None)
  p.add_argument("-m","--multiset_size",
                 help="specify size of multiset",type=int,default=None)
  args = p.parse_args()
  return args

if __name__ == "__main__":
  args = parse_command_line()
  if args.weights_table and args.num_symbols and args.multiset_size:
    weights_table = multiset_weights(args.num_symbols,args.multiset_size)
    print("const unsigned int multiset_weights_table_alphasize = {};"
           .format(args.num_symbols))
    print("const unsigned int multiset_weights_table_width = {};"
           .format(args.multiset_size))
    print("const unsigned int multiset_weights_table[] = {")
    for pos, weights in enumerate(reversed(weights_table)):
      print("  /* pos={} */".format(pos))
      finalsym = "," if pos < len(weights_table) - 1 else "};"
      print("  {}{}".format(',\n  '.join(map(str,weights)),finalsym))
    exit(0)
  if args.verify:
    multisets_verify(args.num_symbols)
  elif args.multiset_size is None:
    for multiset_size in range(2,8+1):
      multiset_output(multiset_size,args.num_symbols,False,False)
  else:
    multiset_output(args.multiset_size,args.num_symbols,args.show,args.perm)
