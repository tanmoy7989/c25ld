import os, sys
import numpy as np
import pickle


# Plot-format
fmt = '.svg'

# --------------------Matplotlib rc---------------------------------------------
import matplotlib
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter, MultipleLocator, MaxNLocator
import matplotlib.pyplot as plt
pad_left = 0.2; pad_bottom = 0.2

matplotlib.rc('lines', linewidth = 3, linestyle = 'solid', 
			  marker = 'None', markersize = 6, markeredgewidth = 0.5)
matplotlib.rc('figure', facecolor = 'w', edgecolor = 'w', figsize = (6.5,6.5))
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rc('xtick', labelsize = 16)
matplotlib.rc('ytick', labelsize = 16)
matplotlib.rc('axes', labelsize = 'large') 	

def prune_axis(axis):
	nbins = 5
	xlim = axis.get_xlim(); ylim = axis.get_ylim()
	delta_x = (xlim[1] - xlim[0])/nbins; delta_y = (ylim[1] - ylim[0])/nbins	
	axis.xaxis.set_major_locator(MultipleLocator(delta_x))
	axis.yaxis.set_major_locator(MultipleLocator(delta_y))

#-------------------------------------------------------------------------------

# data and target locations
c10_dir = os.path.abspath('../data/analysis/aug15_transferability/c10')
c20_dir = os.path.abspath('../data/analysis/aug15_transferability/c20')
c30_dir = os.path.abspath('../data/analysis/aug15_transferability/c30')
c40_dir = os.path.abspath('../data/analysis/aug15_transferability/c40')
c50_dir = os.path.abspath('../data/analysis/aug15_transferability/c50')
#c25_poly_dir = os.path.abspath('../data/analysis/polymer')

savedir = sys.argv[1]
fftype = sys.argv[2]

### ------------------- Plot_dicts----------------------------------------------

#ncases = 6
#case_dirs = {1: c25_poly_dir, 2: c25_dir, 3: c40_dir, 4: c50_dir, 5: c12_dir, 6: c18_dir}
#titles = {1: 'dry_polymer', 2: 'c25', 3: 'c40', 4: 'c50', 5: 'c12 X 3', 6: 'c18 X 2'} 

ncases = 5
case_dirs = {1: c10_dir, 2: c20_dir, 3: c30_dir, 4: c40_dir, 5: c50_dir}
titles = {1: 'c10', 2: 'c20', 3: 'c30', 4: 'c40', 5: 'c50'} 

cgtypes = {1: 'SP', 2: 'SPLD', 3: 'LD'}

nrows = 3 #AA-SP, AA-SPLD, AA-LD
ncols = ncases

axlbls = {'Rg': r'$Rg (\AA)$', 'R_EE' : r'$R_{EE} (\AA)$', 
		  'Dist': 'Distribution', 'pmf' : r'$W(Rg)$' + ' (kJ/mol)',
	      'kappa' : r'$\kappa$', 'SASA' : r'$SASA (\AA^2)$ ' , 
	      'SASA_atom' : r'$\frac{SASA}{atom} (\frac{\AA^2}{atom})$'}

clrs = {1: 'magenta', 2: 'red', 3: 'blue', 4: 'green', 5: 'black', 6: 'purple'}
linestyles = {'AA': '-', 'SP': ':', 'SPLD': '--', 'LD': '-.'}

maketitle = True
def maketitle(axis,j):	
	# use only for visualization
	# do not use while making final figures
	titles = {1: 'Coefficient of shape anisotropy',
			  2: 'Total SASA',
			  3: 'End-to-end distance',
			  4: '',
			  5: 'SASA/atom',
			  6: 'Radius of gyration'
			  }
	
	axis.set_title(titles[j], size = 'xx-large')

#-------------------------------------------------------------------------------


#------------------------ PLotting----------------------------------------------
def plot_hist(ax, pickleName, label, case):
	try:
		data = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'data pickle for %s not found\n' % (pickleName.split('/')[-1].split('.pickle')[0])
		return
	x = data['bin_centers']
	y = data['bin_vals']
	ax.plot(x,y,linestyle = linestyles[label], color = clrs[case], label = label)
	ax.set_xlim([x.min(), x.max()])
	[i.set_color(clrs[case]) for i in ax.get_yticklabels()]


def make1Dplot(geom_prop):
	fig = plt.figure()
	axlist = [0] * (nrows*ncols)
	
	## Plotting
	for case in range(1,ncases+1):
		target_dir = case_dirs[case]
		
		#------------ AA-SP-------------------------------------
		ind = case
		ax = plt.subplot(nrows, ncols, ind)
		axlist[ind-1] = ax
		
		## plot all-atom property
		AA_label = 'AA'
		pickleName = os.path.join(target_dir, '%s_%s_hist1D_%s.pickle' % (AA_label, fftype, geom_prop))
		plot_hist(ax, pickleName, label = AA_label, case = case)
		## plot coarse-grained property
		CG_label = 'SP'
		pickleName = os.path.join(target_dir, 'CG_%s_%s_hist1D_%s.pickle' % (fftype, CG_label, geom_prop))
		plot_hist(ax, pickleName, label = CG_label, case = case)
		
		#------------ AA-SPLD-------------------------------------
		ind = case+ncols
		ax = plt.subplot(nrows, ncols, ind)
		axlist[ind-1] = ax
		
		## plot all-atom property
		AA_label = 'AA'
		pickleName = os.path.join(target_dir, '%s_%s_hist1D_%s.pickle' % (AA_label, fftype, geom_prop))
		plot_hist(ax, pickleName, label = AA_label, case = case)
		## plot coarse-grained property
		CG_label = 'SPLD'
		pickleName = os.path.join(target_dir, 'CG_%s_%s_hist1D_%s.pickle' % (fftype, CG_label, geom_prop))
		plot_hist(ax, pickleName, label = CG_label, case = case)
		
		#------------ AA-LD-------------------------------------
		ind = case+2*ncols
		ax = plt.subplot(nrows, ncols, ind)
		axlist[ind-1] = ax
		
		## plot all-atom property
		AA_label = 'AA'
		pickleName = os.path.join(target_dir, '%s_%s_hist1D_%s.pickle' % (AA_label, fftype, geom_prop))
		plot_hist(ax, pickleName, label = AA_label, case = case)
		## plot coarse-grained property
		CG_label = 'LD'
		pickleName = os.path.join(target_dir, 'CG_%s_%s_hist1D_%s.pickle' % (fftype, CG_label, geom_prop))
		plot_hist(ax, pickleName, label = CG_label, case = case)

	## Axes design
	for i, ax in enumerate(axlist):
		ind = i + 1
		ax.yaxis.set_major_locator(MaxNLocator(nbins = 4, prune = 'both'))
		ax.xaxis.set_major_locator(MaxNLocator(nbins = 4, prune = 'both'))
		if (ind<=ncols*2):
			ax.set_xticklabels([])
		
		if [1,1+ncols,1+2*ncols].__contains__(ind):
			leg = ax.legend(loc = 'best', prop = {'size': 16})
			leg.get_frame().set_alpha(0.4)		
		
		#if (range(2,ncols) + range(2+ncols,2*ncols) + range(2+2*ncols,3*ncols)).__contains__(ind):
		#	ax.set_yticklabels([])
		
	for case in range(1, ncases+1):
		axlist[case-1].set_title(titles[case], color = clrs[case], fontsize = 16, weight = 'bold')
	
	plt.subplots_adjust(wspace = 0.25, hspace = 0., left = pad_left, bottom = pad_bottom)
	
	plt.figtext(0.45, 0.07, axlbls[geom_prop], size = 'xx-large')
	plt.figtext(0.03, 0.45, axlbls['Dist'], size = 'xx-large', rotation = 90)
	
	figname = os.path.join(savedir, geom_prop + '_' + fftype + fmt)
	plt.savefig(figname, dpi = 300)




## MAIN
proplist = ['kappa', 'Rg', 'R_EE', 'SASA', 'SASA_atom']
for prop in proplist:
	make1Dplot(prop)

plt.show()
