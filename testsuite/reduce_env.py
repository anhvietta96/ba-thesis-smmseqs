import argparse
import re

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
                    prog='Qgram Environment Reducer',
                    description='reduce a full environment to a percentage')
  parser.add_argument('filname')
  parser.add_argument('-p','--percentage', type=int,
                    help='an integer for the accumulator')
  
  args = parser.parse_args()
  filename = args.filename

  file1 = open(filename, 'r')
  
  single_qgram_dict = []
  summary_dict = {}
  number_pattern = r'-?\d+'
  sorted_qgram_env = file1.readlines()

  for sorted_qgram_idx, env in enumerate(sorted_qgram_env):
    single_qgram_dict.append({})
    num_list = re.findall(number_pattern,env)
    for idx,num in enumerate(num_list):
      if num not in single_qgram_dict[sorted_qgram_idx].keys():
        single_qgram_dict[sorted_qgram_idx][num] = []
      single_qgram_dict[sorted_qgram_idx][num].append(idx)

      if num not in summary_dict.keys():
        summary_dict[num] = []
      summary_dict.append(idx)

  print(single_qgram_dict)
  print(summary_dict)
  
  file1.close()
