#!/usr/bin/env python3

import sys, argparse, re, os, functools

class Scorematrix:
  def __init__(self,inputfile):
    background_re = re.compile(r'^# Background \(precomputed optional\): ')
    self.scorematrix = dict()
    self.characters = None
    self._maxscore = None
    self._minscore = None
    self._floatmatrix = None
    self._background = None
    self._alp_format = False
    try:
      stream = open(inputfile)
    except IOError as err:
      sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
      exit(1)
    scorelinenum = 0
    for line in stream:
      if re.search(r'^#',line):
        if re.search(background_re,line):
          background_values = re.sub(background_re,'',line)
          self._background = [float(s) for s in background_values.split()]
          total = sum(self._background)
        continue
      if re.search(r'^\s+',line):
        self.characters = re.sub(r'\s+','',line)
        for aa in self.characters:
          self.scorematrix[aa] = dict()
      else:
        l = line.split()
        aa = l[0]
        if aa != self.characters[scorelinenum]:
          sys.stderr.write('l[0] = {} != {} = characters[{}]\n'
                            .format(l[0],self.characters[scorelinenum],
                                    scorelinenum))
          exit(1)
        scorelinenum += 1
        if len(l) != len(self.characters) + 1:
          sys.stderr.write(('{}: line\n{}contains {} scores but must ' +
                            'contain {}\n')
                            .format(sys.argv[0],line,len(l) - 1,
                                    len(self.characters)))
          exit(1)
        for idx in range(1,len(l)):
          if re.search(r'\.',l[idx]):
            assert self._floatmatrix is None or self._floatmatrix
            score = float(l[idx])
            self._floatmatrix = True
          else:
            assert self._floatmatrix is None or not self._floatmatrix
            score = int(l[idx])
            self._floatmatrix = False
          self.update_minmax(score)
          bc = self.characters[idx-1]
          self.scorematrix[aa][bc] = score
  def alp_format_set(self):
    self._alp_format = True
  def __str__(self):
    return self.pretty_print(self.characters)
  def is_float_matrix(self):
    return self._floatmatrix
  def pretty_print(self,characters):
    lines = list()
    if not self._alp_format:
      lines.append('   ' + '  '.join(list(characters)))
    else:
      lines.append('{}'.format(len(characters)))
    for a in characters:
      line = list()
      if not self._alp_format:
        line.append('{}'.format(a))
      for b in characters:
        score = self.scorematrix[a][b]
        if isinstance(score,int):
          if score < 0 or score > 9:
            line.append(' {}'.format(score))
          else:
            line.append('  {}'.format(score))
        else:
          line.append('  {:.4f}'.format(score))
      if self._alp_format:
        this_line = ''.join(line)
        lines.append(re.sub(r'^\s+','',this_line))
      else:
        lines.append(''.join(line))
    return '\n'.join(lines)
  def constexpr(self,name,characters = None):
    if characters is None:
      characters = self.characters
    num_of_chars = len(characters)
    smallest_score = min([self.scorematrix[a][b] for a in characters
                                                 for b in characters])
    score_lines = list()
    for a in characters:
      line = list()
      for b in characters:
        score = self.scorematrix[a][b]
        if score >= 0 and score < 10:
          line.append(' {}'.format(score))
        else:
          line.append('{}'.format(score))
      score_lines.append('/* {} */ {{{}}}'.format(a,','.join(line)))
    pp_var = '{}_HPP'.format(name.upper())
    gen_comment = '/* generated by {} DO NOT EDIT */'.format(' '.join(sys.argv))
    ifdef = '#ifndef {}\n'.format(pp_var)
    define = '#define {}\n'.format(pp_var)
    include = '#include <cstdint>\n'
    decl = 'struct {}{} {{'.format(name[0].upper(),name[1:])
    characters_decl = ('static constexpr const char characters[] = \"{}\";'
                       .format(characters))
    num_of_chars_decl = ('static constexpr size_t num_of_chars = size_t({});'
                          .format(len(characters)))
    character_spec = ('static constexpr const char character_spec[] = \"{}\";'
                      .format('|'.join(list(characters))))
    matrix_decl = ('static constexpr const int8_t score_matrix[{}][{}] = {{'
                   .format(num_of_chars,num_of_chars))
    header_line = '{}/* {} */'.format(' ' * 6,'  '.join(list(characters)))
    smallest = ('static constexpr const int8_t smallest_score = {};'
                .format(smallest_score))
    return gen_comment + ifdef + define + include + decl + characters_decl \
                       + num_of_chars_decl + character_spec + matrix_decl \
                       + header_line + (','.join(score_lines)) \
                       + '};' + smallest + '}; #endif'
  def getscore(self,a,b):
    if a in self.scorematrix and b in self.scorematrix[a]:
      return self.scorematrix[a][b]
    else:
      raise ValueError ('{}: getscore({},{}) is not defined\n'
                        .format(sys.argv[0],a,b))
  def setscore(self,a,b,score):
    assert a in self.scorematrix and b in self.scorematrix[a]
    self.scorematrix[a][b] = score
  def orderscore_get(self,numeric=False):
    def cmp_by_score(p1,p2):
      if p1[1] < p2[1]:
        return 1
      if p1[1] > p2[1]:
        return -1
      return 0
    orderscore = dict()
    for cc_idx, cc in enumerate(self.characters):
      # sort in descending order of score values of the characters replacing
      # cc
      score_line = self.scorematrix[cc]
      assert len(score_line) == len(self.characters)
      this_order_score = list()
      for bb_idx, bb in enumerate(self.characters):
        if numeric:
          this_order_score.append((bb_idx,score_line[bb]))
        else:
          this_order_score.append((bb,score_line[bb]))
      this_order_score = sorted(this_order_score,
                                key=functools.cmp_to_key(cmp_by_score))
      if numeric:
        assert not (cc_idx in orderscore)
        orderscore[cc_idx] = this_order_score
      else:
        assert not (cc in orderscore)
        orderscore[cc] = this_order_score
    return orderscore
  def characters_get(self):
    return self.characters
  def alphasize_get(self):
    return len(self.characters)
  def background_get(self):
    if self._background and len(self._background) != self.alphasize_get() + 1:
      sys.stderr.write(('{}: number of background values {} differs from '
                        'alphabet size {}\n')
                        .format(sys.argv[0],
                                len(self._background),
                                self.alphasize_get()))
      exit(1)
    return self._background
  def eval_score(self,u,v):
    assert len(u) == len(v)
    score = 0
    for a,b in zip(u,v):
      score += self.scorematrix[a][b]
    return score
  def update_minmax(self,score):
    if self._minscore is None or self._minscore > score:
      self._minscore = score
    if self._maxscore is None or self._maxscore < score:
      self._maxscore = score
  def minscore_get(self):
    return self._minscore
  def maxscore_get(self):
    return self._maxscore

def show_verbose(characters,background):
  print('# alphabet_size={}'.format(len(characters)))
  print('# alphabet={}'.format(characters))
  if background:
    background_lines = ['# {} {}'.format(a,b) \
                        for a,b in zip(characters,background)]
    print('\n'.join(background_lines))

def show_alp_background(scorematrixfile,characters,background):
  base = os.path.basename(scorematrixfile)
  outfile = '{}.bcg'.format(base)
  try:
    stream = open(outfile,'w')
  except IOError as err:
    sys.stderr.write('{}: {}\n'.format(sys.argv[0],err))
    exit(1)
  stream.write('{}\n'.format(len(characters)))
  cut_background = background[0:len(characters)]
  total = sum(cut_background)
  addtovalue = (1.0 - total)/float(len(characters))
  total = 0.0
  for f in cut_background:
    new_f = f + addtovalue
    stream.write('{:.10f}\n'.format(new_f))
    total += new_f
  stream.close()

def parse_command_line(argv):
  p = argparse.ArgumentParser()
  p.add_argument('-a','--alphabet',metavar='<string>',
                  help='specify alphabet for which to output the score matrix',
                  type=str,default=None)
  p.add_argument('-v','--verbose',action='store_true',default=False,
                  help='be verbose and provide additional information')
  p.add_argument('--constexpr',type=str,metavar='name',default=None,
                  help=('output matrix in form of constexpr, argument is '
                        'of filename'))
  p.add_argument('--alp',action='store_true',default=False,
                  help=('generate output of matrix and frequency values in '
                        'format compatible for alp'))
  p.add_argument('scorematrix',
                  help='specify file with score_matrix',type=str,
                  default='blosum62-float.txt')
  args = p.parse_args(argv)
  if args.alp and args.verbose:
    sys.stderr.write('{}: option -v and --alp are not compatible\n'
                      .format(sys.argv[0]))
    exit(1)
  return args

if __name__ == '__main__':
  args = parse_command_line(sys.argv[1:])
  scorematrix = Scorematrix(args.scorematrix)
  if args.alp:
    scorematrix.alp_format_set()
  if args.constexpr:
    if args.alphabet:
      print(scorematrix.constexpr(args.constexpr,args.alphabet))
    else:
      print(scorematrix.constexpr(args.constexpr))
  elif args.alphabet:
    if args.verbose:
      show_verbose(args.alphabet,scorematrix.background_get())
    if args.alp:
      show_alp_background(args.scorematrix,
                          args.alphabet,scorematrix.background_get())
      print(scorematrix.pretty_print(args.alphabet))
  else:
    if args.verbose:
      show_verbose(scorematrix.characters_get(),scorematrix.background_get())
    if args.alp:
      show_alp_background(args.scorematrix,
                          scorematrix.characters_get(),
                          scorematrix.background_get())
    print('{}'.format(scorematrix))
