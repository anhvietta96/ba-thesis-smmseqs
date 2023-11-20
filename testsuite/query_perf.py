import matplotlib.pyplot as plt
import numpy as np

groups = ("Sort method without SIMD","Sort method with SIMD","MMseqs2 Method")
measurements = [{
    '1': (679, 695, 3331),
    '10': (675, 705, 3349),
    '100': (1033, 1144, 3715),
    '1000': (3681, 3948, 6210),
},{
    '1': (635, 703, 3432),
    '10': (680, 726, 3391),
    '100': (974, 1014, 3800),
    '1000': (3668, 3840, 6415),
},{
    '1': (667, 661, 3310),
    '10': (718, 730, 3390),
    '100': (1287, 1311, 4154),
    '1000': (5961, 6192, 9019),
}]
ref = 210848
generated_ctxfree = [(162339,1857548,17981137,148986494),(259269,1555134,15615308,139012013),(227180,2420182,24048924,217952274)]
generated_ctxsens = [(162339,1857548,17981137,148986494),(259269,1555134,15615308,139012013),(227180,2420182,24048924,217952274)]
time_generated_ctxsens = [(3642,3672,3679,3671),(3684,4056,3825,3710),(3742,3644,3764,4765)]
max = (7000,7000,10000)
seeds = ('1101011','11101101','1111010101')

#Time to Generate 
'''
for l in range(len(measurements)):
  x = np.arange(len(groups))  # the label locations
  width = 0.2  # the width of the bars
  multiplier = 0

  fig, ax = plt.subplots(layout='constrained')

  for attribute, measurement in measurements[l].items():
      offset = width * multiplier
      rects = ax.bar(x + offset, measurement, width, label=attribute)
      ax.bar_label(rects, padding=3)
      multiplier += 1

  # Add some text for labels, title and custom x-axis tick labels, etc.
  ax.set_ylabel('Time (ms)')
  ax.set_xlabel('Generation method')
  ax.set_title('Generate time of seed ' + seeds[l])
  ax.set_xticks(x + width, groups)
  ax.legend(loc='upper left', ncols=3)
  ax.set_ylim(0, max[l])

  plt.show()
'''

#Estimation evaluation
'''computed = {}
for i in range(len(generated_ctxsens[0])):
  computed[str(10**i)] = []

for seed_measurement in generated_ctxsens:
  for i,measurement in enumerate(seed_measurement):
    computed[str(10**i)].append(measurement/ref/10**i)

x = np.arange(len(seeds))  # the label locations
width = 0.2  # the width of the bars
multiplier = 0

fig, ax = plt.subplots()

for attribute, measurement in computed.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, measurement, width, label=attribute)
    ax.bar_label(rects, padding=3, fmt='%.2f')
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Ratio')
ax.set_xlabel('Seeds')
ax.set_title('Ratio of number of generated q-grams using\n context-sensitive approach vs. expected')
ax.set_xticks(x + width, seeds)
ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 2)

plt.show()
'''

#Local composition bias score correction
'''generated_corr = {
   '0': (1857548,1555134,2420182),
   '0.25': (2790996,2635674,4690658),
   '0.5': (5486423,6250202,12753947),
   '0.75': (11277915,15054689,35612937),
   '1': (23446359,36885976,101611619),
}

x = np.arange(len(seeds))  # the label locations
width = 0.18  # the width of the bars
multiplier = 0

fig, ax = plt.subplots()

for attribute, measurement in generated_corr.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, [measurement[i]/ref/10 for i in range(len(measurement))], width, label=attribute)
    ax.bar_label(rects, padding=3, fmt='%.2f')
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Ratio')
ax.set_xlabel('Seeds')
ax.set_title('Ratio of number of generated q-grams using context-free approach\n with varying degree of local composition bias score correction\n vs. expected')
ax.set_xticks(x + width, seeds)
ax.legend(loc='upper left', ncols=3)
ax.set_ylim(0, 52)

plt.show()
'''

hist = [1.57128*10**-10,
        4.99069*10**-9,
        7.43898*10**-8,
        6.95762*10**-7,
        4.61309*10**-6,
        2.3255*10**-5,
        9.34262*10**-5,
        0.000309625,
        0.000869178,
        0.00211017,
        0.00450491,
        0.0085726,
        0.0147069,
        0.0229693,
        0.0329433,
        0.0437397,
        0.0541808,
        0.0630971,
        0.0696133,
        0.0733071,
        0.0741987,
        0.0726166,
        0.0690323,
        0.0639395,
        0.0578042,
        0.0510595,
        0.0441075,
        0.0373044,
        0.0309343,
        0.0251886,
        0.0201644,
        0.0158826,
        0.0123129,
        0.00939656,
        0.00706093,
        0.00522715,
        0.00381487,
        0.00274653,
        0.00195138,
        0.00136828,
        0.000946764,
        0.000646445,
        0.000435655,
        0.000289904,
        0.00019056,
        0.000123748,
        7.93793*10**-5,
        5.02825*10**-5,
        3.14482*10**-5,
        1.94218*10**-5,
        1.18474*10**-5,
        7.14002*10**-6,
        4.25101*10**-6,
        2.49941*10**-6,
        1.45056*10**-6,
        8.3075*10**-7,
        4.69538*10**-7,
        2.61963*10**-7,
        1.44284*10**-7,
        7.8428*10**-8,
        4.20445*10**-8,
        2.22153*10**-8,
        1.15652*10**-8,
        5.93279*10**-9,
        2.99963*10**-9,
        1.49464*10**-9,
        7.33496*10**-10,
        3.54183*10**-10,
        1.68155*10**-10,
        7.84772*10**-11,
        3.60141*10**-11,
        1.62557*10**-11,
        7.21515*10**-12,
        3.14687*10**-12,
        1.34755*10**-12,
        5.66595*10**-13,
        2.34019*10**-13,
        9.50305*10**-14,
        3.7932*10**-14,
        1.48673*10**-14,
        5.71709*10**-15,
        2.15582*10**-15,
        7.98874*10**-16,
        2.90902*10**-16,
        1.04149*10**-16,
        3.65371*10**-17,
        1.25353*10**-17,
        4.2181*10**-18,
        1.3906*10**-18,
        4.53193*10**-19,
        1.43834*10**-19,
        4.46753*10**-20,
        1.34371*10**-20,
        3.9427*10**-21,
        1.16384*10**-21,
        3.25056*10**-22,
        9.27098*10**-23,
        2.29296*10**-23,
        6.04994*10**-24,
        1.51732*10**-24,
        3.64802*10**-25,
        1.05426*10**-25,
        1.09227*10**-26,
        5.45161*10**-27,
        0,
        3.62396*10**-28]

x = list(range(-28,78))

fig, ax = plt.subplots()

rects = ax.bar(x, hist[:28] + [0]*(78+28-28),color='lightblue')
rects2 = ax.bar(x, [0]*28 + hist[28:],color='lightyellow')
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Frequency')
ax.set_xlabel('Score')
ax.set_title('Score distribution of heptagram vs heptagram,\nboth uniformly distributed')
ax.set_ylim(0, 0.1)

plt.show()