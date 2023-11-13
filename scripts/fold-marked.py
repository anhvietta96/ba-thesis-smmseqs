#!/usr/bin/env python3

import sys, re, argparse, shutil

def parse_arguments():
  p = argparse.ArgumentParser(description=('fold file in regions marked by '
                                           '%%%BEGIN{fold} END{fold}'))
  p.add_argument('-w','--width',type=int,default=75,
                 help='specify width of folded lines, default 75')
  p.add_argument('--inplace', action='store_true',
                  help='overwrite input file')
  p.add_argument('inputfile',type=str,
                  help='specify input file')
  return p.parse_args()

def print_folded(line,width, output_stream):
  words = re.findall(r'\S+',line)
  current_length = 0
  out_words = list()
  for word in words:
    if current_length + len(out_words) + len(word) > width:
      print(' '.join(out_words), file=output_stream)
      current_length = 0
      out_words.clear()
    out_words.append(word)
    current_length += len(word)
  if out_words:
    print(' '.join(out_words), file=output_stream)

def fold_file(filename,width, inplace):
  if inplace:
    output_stream = open('{}.folded'.format(filename),'w')
  else:
    output_stream = sys.stdout
  try:
    stream = open(filename)
  except IOError as err:
    sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
    exit(1)
  dofold = False
  first_line = True
  newline = False
  for line_num, line in enumerate(stream,1):
    if re.search(r'^%%%END\{fold\}',line):
      if not dofold:
        sys.stderr.write('{}: line {}: END has no BEGIN\n'
                          .format(sys.argv[0],line_num))
        exit(1)
      dofold = False
      print('', file=output_stream)
      newline = False
    elif re.search(r'^%%%BEGIN\{fold\}',line):
      if first_line:
        first_line = False
      else:
        if dofold:
          sys.stderr.write('{}: line {}: missing previous END\n'
                          .format(sys.argv[0],line_num))
          exit(1)
      dofold = True
      if not newline:
        print('', file=output_stream)
    elif dofold:
      print_folded(line,width, output_stream)
      if not re.search('^$',line.rstrip()):
        newline = True
    else:
      print(line,end='', file=output_stream)
  stream.close()
  if inplace:
    output_stream.close()
    shutil.move('{}.folded'.format(filename), filename)

args = parse_arguments()
fold_file(args.inputfile,args.width, args.inplace)
