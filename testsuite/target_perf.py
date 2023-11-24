import matplotlib.pyplot as plt
import numpy as np

'''groups = ("11101","111011","1101011","110010000101","11101101","1101010011","1111010101","11010110011")
measurements = [{
    'Linear hashing': (324, 327, 324, 330, 330, 334, 334, 333),
    'Recursive hashing': (346, 348, 369, 396, 369, 393, 390, 393)
}]
ref = 210848
generated_ctxfree = [(162339,1857548,17981137,148986494),(259269,1555134,15615308,139012013),(227180,2420182,24048924,217952274)]
generated_ctxsens = [(162339,1857548,17981137,148986494),(259269,1555134,15615308,139012013),(227180,2420182,24048924,217952274)]
time_generated_ctxsens = [(3642,3672,3679,3671),(3684,4056,3825,3710),(3742,3644,3764,4765)]
max = (7000,7000,10000)
seeds = ('1101011','11101101','1111010101')'''

#Time to Generate
#for l in range(len(measurements)):
'''x = np.arange(len(groups))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in measurements[0].items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Time (ms)')
ax.set_xlabel('Seed')
#ax.set_title('Generate time of seed ' + seeds[l])
ax.set_xticks(x + width/2, groups)
ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 450)

plt.show()'''
'''
groups = ('sorted-unsorted score matrix','unsorted-unsorted score matrix')
measurements = {
  'q=2': (3,6),
  'q=3': (586,3264)
}
x = np.arange(len(groups))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in measurements.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Time (ms)')
ax.set_xlabel('Matrix type')
#ax.set_title('Generate time of seed ' + seeds[l])
ax.set_xticks(x + width/2, groups)
ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 3500)

plt.show()
'''

'''
groups = ('1','2.5','5','7.5')
measurements = {
  'generated': (1.94,5.792,28.722,231.11)
}
x = np.arange(len(groups))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in measurements.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('q-grams per position / sensitivity')
ax.set_xlabel('Sensitivity')
#ax.set_title('Generate time of seed ' + seeds[l])
ax.set_xticks(x, groups)
#ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 250)

plt.show()
'''

'''
groups = ('1101011','11101101','1111010101')
measurements = {
  'Context-free': (1,1,1),
  'Context-sensitive': (406,475,406)
}
x = np.arange(len(groups))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in measurements.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Time (ms)')
ax.set_xlabel('Seeds')
#ax.set_title('Generate time of seed ' + seeds[l])
ax.set_xticks(x + width/2, groups)
ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 800)

plt.show()
'''

groups = ('10','100','1000')
measurements = [{
  'Target processing': (153,153,153),
  'Query processing': (88,463,4620),
  'Merging data': (20,84,384)
},{
  'Target processing': (175,175,175),
  'Query processing': (140,832,7933),
  'Merging data': (46,79,305)
}]
x = np.arange(len(groups))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0
ylim = [5300,9000]

fig, ax = plt.subplots(layout='constrained')

for i in range(2):
  for attribute, measurement in measurements[1].items():
      offset = width * multiplier
      rects = ax.bar(x + offset, measurement, width, label=attribute)
      ax.bar_label(rects, padding=3)
      multiplier += 1

  # Add some text for labels, title and custom x-axis tick labels, etc.
  ax.set_ylabel('Time (ms)')
  ax.set_xlabel('Sensitivity')
  #ax.set_title('Generate time of seed ' + seeds[l])
  ax.set_xticks(x + width/2, groups)
  ax.legend(loc='upper left', ncols=3)
  ax.set_ylim(0, ylim[1])

  plt.show()
