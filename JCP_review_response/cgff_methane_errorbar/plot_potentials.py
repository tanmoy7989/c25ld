#!/usr/bin/env python

import numpy as np
import parse_potential as pp
import matplotlib.pyplot as plt

runs = [1,2,3,5]
clrs = ['red', 'blue', 'green', 'cyan', 'black']
fig = plt.figure(figsize = (4,4), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)

r_set = []
p_set = []

for run in runs:
    sumfile = 'm25_%d_SPLD_sum.txt' % run
    r, p = pp.parseLog(sumfile)['LD']
    p -= np.max(p)
    r_set.append(r)
    p_set.append(p)

r = np.array(r_set)
p = np.array(p_set)
mu = np.mean(p, axis = 0)
sigma = np.std(p, axis = 0, ddof = 1)

ax.errorbar(r[0,:], mu, yerr = sigma)
ax.set_xlim([0, 24])
ax.set_ylim([-1.2, 0.1])
ax.set_xlabel('local density of methane', fontsize = 10)
ax.set_ylabel('potential (kcal/mol)', fontsize = 10)
ax.legend(loc = 'best')

plt.show()
