#!/usr/bin/env python3

import sys, argparse

def parse_command_line(argv):
  p = argparse.ArgumentParser(description='compare two sets of hits')
  p.add_argument('-d','--debug',action='store_true',default=False,
                  help='show debug output')
  p.add_argument('--with_bit_score',action='store_true',default=False,
                  help='track scoress of hits if available')
  p.add_argument('--half_missing',action='store_true',default=False,
                  help='show half missing values')
  p.add_argument('--pairs',type=int,default=20000 * 426,
                  help=('specify number of possible pair to obatin filter '
                        'efficiency'))
  p.add_argument('gold_standard',type=str,
                  help='specify input file with gold standard')
  p.add_argument('prediction',type=str,
                  help='specify input file with prediction')
  return p.parse_args(argv)

def read_hits(debug,with_bit_score,filename):
  try:
    stream = open(filename)
  except IOError as err:
    sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
    exit(1)
  ids2bit_scores = dict()
  for line in stream:
    line = line.rstrip()
    if line.startswith('# TIME'):
      print(line)
      continue
    elif line.startswith('#'):
      continue
    arr = line.split('\t')
    ids = (arr[0],arr[1])
    if ids not in ids2bit_scores:
      ids2bit_scores[ids] = list()
    if with_bit_score:
      ids2bit_scores[ids].append(float(arr[-1]))
    else:
      ids2bit_scores[ids].append(0)
  stream.close()
  if debug:
    for ids, bit_scores in ids2bit_scores.items():
      print('{}\t{}'.format(ids,bit_scores))
  return ids2bit_scores

def cmp_hits(num_pairs,track_half_missing,gold_standard,prediction):
  for p_ids, p_bit_scores in prediction.items():
    if p_ids not in gold_standard and max(p_bit_scores) >= 56:
      sys.stderr.write('{}: no hit for p_ids {} with {} in gold_standard\n'
                       .format(sys.argv[0],p_ids,p_bit_scores))
      exit(1)
    if p_ids in gold_standard:
      g_bit_scores = gold_standard[p_ids]
      assert len(g_bit_scores) == 1, \
             'g_bit_scores[{}]={}'.format(p_ids,g_bit_scores)
      for bs in p_bit_scores:
        if bs > g_bit_scores[0] + 1:
          diff = 100.0 * (g_bit_scores[0] + 1 - bs)/(g_bit_scores[0] + 1)
          if diff > 2.0:
            sys.stderr.write(('{}: hit score {} for p_ids {} is more than 10% '
                              'larger than gold_standard bit score {}\n')
                               .format(sys.argv[0],bs,p_ids,g_bit_scores[0]))
            exit(1)
  true_positives = 0
  half_missing = dict()
  for g_ids, g_bit_scores in gold_standard.items():
    if g_ids in prediction:
      true_positives += 1
      if track_half_missing:
        bs_max = max(prediction[g_ids])
        diff = int(round(g_bit_scores[0] + 1 - bs_max))
        if diff not in half_missing:
          half_missing[diff] = 0
        half_missing[diff] += 1
  gl = len(gold_standard)
  print('pairs_in_gold_standard={}'.format(gl))
  print('pairs in prediction={}'.format(len(prediction)))
  print('true_positives (predictions present in gold_standard)={}'
        .format(true_positives))
  print('false_negatives (i.e missing predictions)={}'
         .format(gl-true_positives))
  print('sensitivity={:.2f}'.format(100 * (true_positives/gl)))
  print('specificity={:.2f}'.format(100 * (true_positives/len(prediction))))
  print('filter_efficiency (percent)={:.2f}'
        .format(100.0 * len(prediction)/num_pairs))
  if args.half_missing:
    print('half_missing:')
    all_half = sum(half_missing.values())
    for k in sorted(half_missing,key=lambda x:half_missing[x],reverse=True):
      print('{}\t{}\t{:.2f}'.format(k,half_missing[k],
                                    100.0 * half_missing[k]/all_half))

args = parse_command_line(sys.argv[1:])
print('gold_standard')
gold_standard = read_hits(args.debug,args.with_bit_score,args.gold_standard)
print('prediction')
prediction = read_hits(args.debug,args.with_bit_score,args.prediction)
cmp_hits(args.pairs,args.half_missing,gold_standard,prediction)
