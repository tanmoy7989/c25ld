import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#plot format
fmt = '.svg'


#single axes plot
def plot_hist(ax, pickleName, trajtype, label = ''):
	try:
		data = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'data pickle for %s not found\n' % (pickleName.split('/')[-1].split('.pickle')[0])
		return
	x = data[0]
	y = data[1]
	ax.plot(x,y,linestyle = styles[trajtype][0], marker = styles[trajtype][1],
	        markersize = 5, linewidth = 3, color = clrs[trajtype], label = label)
	


#data and target locations
alt_suffix = ''
geom_prop_dir = '/home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop'
c25_dir = os.path.join(geom_prop_dir, 'c25%s/hist' % alt_suffix) #reference
c10_dir = os.path.join(geom_prop_dir, 'c10%s/hist' % alt_suffix)
c20_dir = os.path.join(geom_prop_dir, 'c20%s/hist' % alt_suffix)
c30_dir = os.path.join(geom_prop_dir, 'c30%s/hist' % alt_suffix)
c40_dir = os.path.join(geom_prop_dir, 'c40%s/hist' % alt_suffix)
c50_dir = os.path.join(geom_prop_dir, 'c50%s/hist' % alt_suffix)
fftype = sys.argv[1]
savedir = sys.argv[2]


#plot 
dirs = {0: c25_dir, 1: c10_dir, 2: c20_dir, 3: c30_dir, 4: c40_dir, 5: c50_dir}
clen = {0: 'c25', 1: 'c10', 2: 'c20', 3: 'c30', 4: 'c40', 5: 'c50'} 
trajtypes = ['AA', 'SP', 'SPLD', 'LD']
nclen = len(clen)


#labels colors and styles
axlbls = {'Rg': r'$Rg (\AA)$', 'R_EE' : r'$R_{EE} (\AA)$', 
		  'Dist': 'Distribution', 'pmf' : r'$W(Rg)$' + ' (kJ/mol)',
	      'kappa' : r'$\kappa$', 'SASA_atom' : r'$\frac{SASA}{atom} (\frac{\AA^2}{atom})$'}

clrs = {'AA': 'red', 'SP': 'blue', 'SPLD': 'blue', 'LD': 'blue'}
styles = {'AA': ('None', 'o'), 
          'SP': ('-', 'None'), 
          'SPLD': ('-','None'),
          'LD': ('-', 'None')}


#plot all of 'em together
axlist = [0] * nclen * 3

def make1Dplot(geom_prop):
	fig = plt.figure(figsize = (12,6), facecolor = 'w', edgecolor = 'w')
	currax = 1
	for cgtype in trajtypes[1:]:
	   for n in range(nclen):
	       ax = fig.add_subplot(3,nclen,currax)
	       axlist[currax-1] = ax
	       
	       AA_pickle = os.path.join(dirs[n], '%s_%s_AA.pickle' % (geom_prop, fftype))
	       CG_pickle = os.path.join(dirs[n], '%s_%s_%s.pickle' % (geom_prop, fftype, cgtype))
	       
	       if currax <=6: title = clen[n]
	       else: title = ''
	       if currax == 1: label = 'AA-SP'
	       if currax == 7: label = 'AA-SPLD'
	       if currax == 13: label = 'AA-LD'
	      
	       plot_hist(ax, AA_pickle, trajtype = 'AA')
	       plot_hist(ax, CG_pickle, trajtype = cgtype, label = label)
	       ax.set_title(title)
	       
	       currax += 1
	      
	#design
	currax = 0
	for n, ax in enumerate(axlist):
	   ax.xaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
	   ax.yaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
	   currax = n+1
	   if currax==1 or currax==7 or currax==13: 
	       leg = ax.legend(loc = 'best', prop = {'size': 12})
	       leg.get_frame().set_alpha(0.2)
	   	
	plt.subplots_adjust(hspace = 0.2, wspace = 0.2)      
	figname = os.path.join(savedir, geom_prop + '_' + fftype + fmt)
	plt.savefig(figname, dpi = 300)


#MAIN
proplist = ['kappa', 'Rg', 'R_EE', 'SASA_atom']
for prop in proplist:
	make1Dplot(prop)

plt.show()
