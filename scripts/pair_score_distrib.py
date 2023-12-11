#!/usr/bin/env python3

from Blosum62 import Blosum62Matrix
from math import sqrt

def combine_score_dists(score_dist1,p_score_dist):
  score_dist = dict()
  for p_score, counts in p_score_dist.items():
    for score1, counts1 in score_dist1.items():
      this_score = p_score + score1
      if this_score not in score_dist:
        score_dist[this_score] = 0
      score_dist[this_score] += counts * counts1
  return score_dist

def score_dist1_get(scorematrix):
  score_dist1 = dict()
  for a in scorematrix.alpha:
    for b in scorematrix.alpha:
      this_score = scorematrix.getscore(a,b)
      if this_score not in score_dist1:
        score_dist1[this_score] = 0
      score_dist1[this_score] += 1
  return score_dist1

def show_score_dist(score_dist):
  sum_count = 0
  for s in sorted(score_dist):
    print('{}\t{}'.format(s,score_dist[s]))
    sum_count += score_dist[s]
  alpha_size = int(sqrt(sum_count))
  assert alpha_size * alpha_size == sum_count
  print('sum_count\t{}'.format(sum_count))

bm = Blosum62Matrix()
score_dist1 = score_dist1_get(bm)
print('score_dist1')
show_score_dist(score_dist1)
score_dist_list = [score_dist1]
for q in range(2,7+1):
  score_dist = combine_score_dists(score_dist1,score_dist_list[q-2])
  print('score_dist{}'.format(q))
  show_score_dist(score_dist)
  score_dist_list.append(score_dist)
