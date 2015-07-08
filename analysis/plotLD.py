import numpy as np
import matplotlib.pyplot as plt
import pickle

fig = plt.figure(figsize = (4,4), facecolor = 'w', edgecolor = 'w')
ax = fig.add_subplot(1,1,1)

(AA, SP, SPLD) = pickle.load(open('LD_hist.pickle', 'r'))

lbls = ['AA', 'SP', 'SPLD']
clrs = ['red', 'blue', 'green']
linestyles = ['None', 'solid', 'solid']
markers = ['o', 'None', 'None']
for i, h in enumerate([AA,SP, SPLD]):	
	r = h[0]
	hist = h[1]
	ax.plot(r,hist,linestyle = linestyles[i], marker = markers[i], linewidth = 3, markersize = 6,
		color = clrs[i], label = lbls[i])

ax.legend(loc = 'best', prop = {'size' : 15}, fancybox = True, shadow = True)
ax.set_xlabel('local density of monomers', fontsize = 'large')
ax.set_ylabel('distribution', fontsize = 'large')
plt.subplots_adjust(left = 0.2, bottom = 0.2)
plt.show() 
	
