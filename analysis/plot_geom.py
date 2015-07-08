import os, sys
import numpy as np
import pickle

# --------------------- Plot-dicts----------------------------------------------
axlbls = {'Rg': r'$Rg (\AA)$', 'R_EE' : r'$R_{EE} (\AA)$', 
		  'Dist': 'Distribution', 'pmf' : r'$W(Rg)$' + ' (kJ/mol)',
	      'kappa' : r'$\kappa$', 'SASA' : r'$SASA (\AA^2)$ ' , 
	      'SASA_atom' : r'$SASA/atom (\AA^2/atom)$'}

clrs = {'AA_lj': 'red', 'CG_lj_SP' : 'blue', 'CG_lj_SPLD' : 'green',
		'CG_lj_LD': 'cyan',
		'AA_wca': 'red', 'CG_wca_SP': 'blue', 'CG_wca_SPLD': 'green', 
		'CG_wca_LD': 'cyan'}

lbls_lj = ['AA_lj', 'CG_lj_SP', 'CG_lj_SPLD', 'CG_lj_LD']
lbls_wca = ['AA_wca', 'CG_wca_SP', 'CG_wca_SPLD', 'CG_wca_LD']

#-------------------------------------------------------------------------------

# Plot-format
fmt = '.svg'

# --------------------Matplotlib rc---------------------------------------------
import matplotlib
import matplotlib.cm as cm
from matplotlib.ticker import NullFormatter, MultipleLocator, MaxNLocator
import matplotlib.pyplot as plt
pad_left = 0.2; pad_bottom = 0.2

matplotlib.rc('lines', linewidth = 2, linestyle = 'solid', 
			  marker = 'None', markersize = 6, markeredgewidth = 0.5)
matplotlib.rc('figure', facecolor = 'w', edgecolor = 'w', figsize = (5,5))
matplotlib.rcParams['xtick.major.size'] = 6
matplotlib.rcParams['ytick.major.size'] = 6
matplotlib.rc('xtick', labelsize = 18)
matplotlib.rc('ytick', labelsize = 18)
matplotlib.rc('axes', labelsize = 'xx-large') 	

def prune_axis(axis):
	nbins = 5
	xlim = axis.get_xlim(); ylim = axis.get_ylim()
	delta_x = (xlim[1] - xlim[0])/nbins; delta_y = (ylim[1] - ylim[0])/nbins	
	axis.xaxis.set_major_locator(MultipleLocator(delta_x))
	axis.yaxis.set_major_locator(MultipleLocator(delta_y))

maketitle_flag = False
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
#--------------------------------------------------------------------------------


##### conversion factor
kB = 1.3806e-23
Temp = 298.0
N_A = 6.023e23
conv_factor = (kB*Temp) * N_A/1000.0  # converting from kT to kJoule/mol

# data and target locations
AA_dir = sys.argv[1]
CG_dir = sys.argv[2]
savedir = sys.argv[3]
fftype = sys.argv[4]

if not CG_dir:	CG_dir = AA_dir
if fftype == 'lj':	lbls = lbls_lj
else:	lbls = lbls_wca

def get_target_dir(trjtype):
	if ['AA_lj', 'AA_wca'].__contains__(trjtype):
		target_dir = AA_dir
	elif ['CG_lj_SP', 'CG_lj_SPLD', 'CG_lj_LD', 'CG_wca_SP', 'CG_wca_SPLD', 'CG_wca_LD'].__contains__(trjtype):
		target_dir = CG_dir 
	
	return target_dir


### ------------------- Plotting------------------------------------------------

# 1) Distributions of kappa
fig = plt.figure()
ax =  plt.subplot('111')
for trjtype in lbls:
	target_dir = get_target_dir(trjtype)
	pickleName = os.path.join(os.path.abspath(target_dir), trjtype + '_hist1D_kappa.pickle')
	try:
		hist = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'kappa Pickle for for %s not found' % trjtype
		continue
		
	x = hist['bin_centers']
	y = hist['bin_vals']
	ax.plot(x,y, linestyle = 'solid', color = clrs[trjtype], label = trjtype)

ax.set_xlabel(axlbls['kappa']); ax.set_ylabel(axlbls['Dist'])
ax.legend(loc = 'best', shadow = True)
prune_axis(ax)
if maketitle_flag: maketitle(ax,1)
figname = os.path.join(savedir, 'kappa_'+ trjtype.split('_')[1]+ fmt)
plt.subplots_adjust(left = pad_left, bottom = pad_bottom)
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)


# 2) Distribution of total SASA
fig = plt.figure()
ax =  plt.subplot('111')
for trjtype in lbls:
	target_dir = get_target_dir(trjtype)
	pickleName = os.path.join(os.path.abspath(target_dir), trjtype + '_hist1D_SASA.pickle')
	try:	
		hist = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'SASA pickle for %s not found' % trjtype
		continue
		
	x = hist['bin_centers']
	y = hist['bin_vals']
	ax.plot(x,y, linestyle = 'solid', color = clrs[trjtype], label = trjtype)
	
ax.set_xlabel(axlbls['SASA']); ax.set_ylabel(axlbls['Dist'])
ax.legend(loc = 'best', shadow = True)
prune_axis(ax)
if maketitle_flag: maketitle(ax,2)
figname = os.path.join(savedir, 'SASA_'+ trjtype.split('_')[1]+ fmt)
plt.subplots_adjust(left = pad_left, bottom = pad_bottom)
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)



# 3) Distribution of R_EE
fig = plt.figure()
ax =  plt.subplot('111')
for trjtype in lbls:
	target_dir = get_target_dir(trjtype)
	pickleName = os.path.join(os.path.abspath(target_dir), trjtype + '_hist1D_R_EE.pickle')
	try:	
		hist = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'R_EE pickle for %s not found' % trjtype
		continue
		
	x = hist['bin_centers']
	y = hist['bin_vals']
	plt.plot(x,y, linestyle = 'solid', color = clrs[trjtype], label = trjtype)

ax.set_xlabel(axlbls['R_EE']); ax.set_ylabel(axlbls['Dist'])
ax.legend(loc = 'best', shadow = 'True')
prune_axis(ax)
if maketitle_flag: maketitle(ax,3)
figname = os.path.join(savedir, 'R_EE_'+ trjtype.split('_')[1]+ fmt)
plt.subplots_adjust(left = pad_left, bottom = pad_bottom)
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)



# 4) 2D PMF of Rg vs R_EE
fig = plt.figure()
index = '22'
for i, trjtype in enumerate(lbls):
	target_dir = get_target_dir(trjtype)
	pickleName = os.path.join(os.path.abspath(target_dir), trjtype + '_hist2D_Rg_R_EE.pickle')
	try:
		hist = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'Rg_R_EE 2D pickle for %s not found' %  trjtype
		continue
		
	ax = plt.subplot(index+str(i+1))
	Rg = hist['bin_centers'][0]
	R_EE = hist['bin_centers'][1]
	pmf = hist['pmf']
	X = (Rg - Rg.min())/(Rg.max() - Rg.min())
	Y = (R_EE - R_EE.min())/(R_EE.max() - R_EE.min())
	Z = pmf * conv_factor	
	im = ax.imshow(Z, cmap = cm.spectral_r, 
				   extent = [X.min(),X.max(),Y.min(),Y.max()])
				   	
	## axis labels and annotation
	annotate_txt = trjtype.split('_')[-1]
	if annotate_txt == 'lj' or annotate_txt == 'wca':
		annotate_txt = 'AA'
	ax.annotate(annotate_txt, xy = (0.2, 0.2), fontsize = 18, weight = 'bold')
	ax.xaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
	ax.yaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
	ax.set_xlabel(''); ax.set_ylabel('')
	if i == 0 or i == 1: 
		ax.set_xticklabels([])
	if i == 1 or i == 3:
		ax.set_yticklabels([])

# colorbar
cax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
cbar = fig.colorbar(im, cax = cax, orientation = 'vertical')
cbar.set_label(r'$pmf(Rg, R_{EE}) (kJ/mol)$', size = 'large', rotation = 90)
# central x and y labels 	
plt.figtext(0.35, 0.03, 'Dimensionless ' + r'$Rg$', size = 'x-large')
plt.figtext(0.03, 0.35, 'Dimensionless ' + r'$R_{EE}$', size = 'x-large', rotation = 90)

plt.subplots_adjust(wspace = 0., hspace = 0., left = pad_left, bottom = pad_bottom) #, right = 1 - pad_left)   
figname = os.path.join(savedir, 'Rg_vs_R_EE_'+ trjtype.split('_')[1]+ fmt)
plt.savefig(figname, dpi = 300)


# 5) Distribution of per atom SASA
fig = plt.figure()
ax =  plt.subplot('111')
for trjtype in lbls:
	target_dir = get_target_dir(trjtype)
	pickleName = os.path.join(os.path.abspath(target_dir), trjtype + '_hist1D_SASA_atom.pickle')
	try:
		hist = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'SASA_atom pickle for %s not found' % trjtype
		continue
		
	x = hist['bin_centers'][2:-2]
	y = hist['bin_vals'][2:-2]
	ax.plot(x,y, linestyle = 'solid', color = clrs[trjtype], label = trjtype)
	
ax.set_xlabel(axlbls['SASA_atom']); ax.set_ylabel(axlbls['Dist'])
ax.legend(loc = 'best', shadow = True)
figname = os.path.join(savedir, 'SASA_atom_'+ trjtype.split('_')[1]+ fmt)
prune_axis(ax)
if maketitle_flag: maketitle(ax,5)
plt.subplots_adjust(left = pad_left, bottom = pad_bottom)
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)


# 6) Distributions of Rg from unbiased samples
fig = plt.figure()
ax = plt.subplot('111')
for trjtype in lbls:
	target_dir = get_target_dir(trjtype)
	pickleName = os.path.join(os.path.abspath(target_dir), trjtype + '_hist1D_Rg.pickle')
	try:
		hist = pickle.load(open(pickleName, 'r'))
	except IOError:
		print 'Rg pickle for %s not found' % trjtype
		continue
		
	x = hist['bin_centers']
	y = hist['bin_vals']
	ax.plot(x,y, linestyle = 'solid', color = clrs[trjtype], label = trjtype)
	
ax.set_xlabel(axlbls['Rg']); ax.set_ylabel(axlbls['Dist'])
ax.legend(loc = 'best', shadow = True)
prune_axis(ax)
if maketitle_flag: maketitle(ax,6)
figname = os.path.join(savedir, 'Rg_'+ trjtype.split('_')[1]+ fmt)
plt.subplots_adjust(left = pad_left, bottom = pad_bottom)
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)


plt.show()
