import os, pickle, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import parse_potential as pp
'''
This script makes a set of plots
due for the nsf report submitted on 23rd March, 2015
'''

savedir = 'nsf_Report_plots'
if not os.path.isdir(savedir):
	os.mkdir(savedir)

kB = 1.3806e-23
Temp = 298.0
N_A = 6.023e23
conv_factor = (kB*Temp) * N_A/1000.0  # converting from kT to kJoule/mol

# set matplotlib rc for better quality plots
def setup():
	mpl.rc('xtick', labelsize = 16)
	mpl.rc('ytick', labelsize = 16)
	mpl.rcParams['xtick.major.pad'] = 4
	mpl.rcParams['ytick.major.pad'] = 4
	mpl.rc('font', **{'size': 16, 'weight': 'medium'})
	mpl.rc('lines', linewidth = 3, markersize = 8)
	mpl.rc('figure', facecolor = 'white', edgecolor = 'white')
	mpl.rc('legend', fancybox = True, shadow = True, fontsize = 15)

	
def plot1():
	'''
	Plot sensitivity of Local Density to
	first shell waters for c25ld
	
	size : height = 2.7 inch, width = 7.9 inch
	'''
	setup()
	fig = plt.figure(figsize = (8.2, 3.2))
	plt.subplots_adjust(bottom = 0.22, wspace = 0.005)
	
	ax1 = plt.subplot(121)
	datafile = os.path.abspath('./data/analysis/feb15_runs/rdf/firstshell_lj.pickle')
	data = pickle.load(open(datafile, 'r'))
	ind = data[0].index(6.5)
	nw = data[2]
	ld = data[1][:,ind]
	ax1.plot(ld, nw, linestyle = 'none', marker = 'x', markersize = 3.5, 
			markerfacecolor = 'red', markeredgecolor = 'red',
			label = 'cut = ' + r'$6.5 \AA$')
	ax1.set_xlabel('local density of monomer')
	ax1.set_ylabel('first shell waters')
	ax1.set_xticklabels(['', 5, 10, 15, 20])
	ax1.legend(loc = 'best')

	
	axL2 = plt.subplot(122)
	ax2 = axL2.twinx()
	datafile = os.path.abspath('./data/analysis/feb15_runs/rdf/corrcoeff.pickle')
	data = pickle.load(open(datafile, 'r'))
	ldcuts = data[0]
	corrcoeffsq = data[1][:,0]**2.
	ax2.plot(ldcuts, corrcoeffsq, linestyle = 'solid', linewidth = 3, color = 'blue',
			marker = 'o', markerfacecolor = 'red', markeredgecolor = 'black')
	ax2.plot(6.5*np.ones(10), np.linspace(0, ax2.get_ylim()[1], 10),
			linestyle = 'dashed', linewidth = 2, color = 'black')
	axL2.set_xlabel('local density cut ' + r'$(\AA)$')
	axL2.set_xticklabels(['', 2,4,6,8,10])
	axL2.set_yticklabels('')
	ax2.set_ylabel(r'$R^2$')
	ax2.set_yticklabels([0, '', 0.2, '', 0.4, '', 0.6, '', 0.8])
	ax2.annotate('cut = ' + r'$6.5 \AA$', xy = (6.8, 0.05))
	figname = os.path.join(savedir, 'plot1.png') 
	plt.savefig(figname, bbox_inches = 'tight', pad_inches = 0., dpi = 300)

	
	

def plot2():
	'''
	Plot the LD potential and LD distribution without spline
	for LJ case
	'''	
	setup()
	figsize = (4, 3.4)
	fig = plt.figure(figsize = figsize)
	plt.subplots_adjust(left = 0.22, bottom = 0.22, right = 0.78)
	
	ax = plt.subplot(111)
	target_dir = os.path.abspath('./data/cg_ff/feb15_runs_fsw')	
	prefix = os.path.join(target_dir, 'c25_lj_LD')
	logfile, histfile = pp.parseFileNames(prefix)
	log_data = pp.parseLog(logfile)['LD']
	ld = log_data[0]
	pld = log_data[1]
	hist_data = pp.parseHist(histfile, '', [])['LD']
	ld_t = hist_data[0]; h_t = hist_data[1]; ld_m = hist_data[2]; h_m = hist_data[3]
	ax.plot(ld, pld, linestyle = 'solid', linewidth = 3, color = 'red')
	
	ax2 = ax.twinx()
	ax2.plot(ld_t, h_t, linestyle = 'solid', linewidth = 3, color = 'black', label = 'AA')
	ax2.plot(ld_m,h_m, linestyle = 'dashed', linewidth = 3, color = 'blue', label = 'CG')
	
	ax.set_xlabel('Local Density', fontsize = 13)
	ax.set_ylabel('LD Potential', fontsize = 13, color = 'red')
	ax2.set_ylabel('LD Distribution', fontsize = 13, color = 'blue')
	ax.set_yticks([-0.2, -0.1, 0.0, 0.1, 0.2])
	ax2.set_yticks([0.02, 0.04, 0.06, 0.08, 0.10])
	
	leg2 = ax2.legend(loc = 'best', prop = {'size': 10}, fancybox = False, shadow = False)
	leg2.get_frame().set_alpha(0.4)
	
	figname = os.path.join(savedir, 'plot2.svg') 
	plt.savefig(figname, bbox_inches = 'tight', dpi = 300, pad_inches = 0)


def plot3():
	'''
	plot the LD distribution with and without
	LD potentials for TIP4P water
	'''
	rhos = [1020, 1100, 1170]; Temp = 300
	clrs = ['red', 'blue', 'green']
	ffdir = os.path.expanduser('~/TIP4P/cg_ff')
	
	setup()
	tickloc = mpl.ticker.MaxNLocator(nbins = 8, prune ='both')
	fig = plt.figure(figsize = (4.55,3.7))
	plt.subplots_adjust(left = 0.13, bottom = 0.11, hspace = 0.005)
	
	ax1 = fig.add_subplot(211)
	cgtype = 'SP'
	for i, rho in enumerate(rhos):
		prefix = os.path.join(ffdir, 'R%dT%d_%s' % (rho, Temp,cgtype))
		log, hist = pp.parseFileNames(prefix)
		hist = pp.parseHist(hist, '', []) 
		hist = hist['LD']
		r_t = hist[0]; h_t = hist[1]; r_m = hist[2]; h_m = hist[3]
		ax1.plot(r_t+i, h_t, linestyle = 'solid', color = clrs[i], 
				label = r'$\rho = $' + str(rho) + ' '+ r'$kg/m^3$')
		ax1.plot(r_m+i, h_m, linestyle = 'dashed', color = clrs[i])
		ax1.set_ylabel('LD Distribution')
		leg1 = ax1.legend(loc = 2)
		leg1.get_frame().set_alpha(0.05)
		ax1.annotate('Pair spline\nonly', xy = (18,0.20))
		ax1.set_xticklabels('', visible = False)
		ax1.yaxis.set_major_locator(tickloc)
	
	ax2 = fig.add_subplot(212, sharex = ax1)
	cgtype = 'SPLD'
	for i, rho in enumerate(rhos):
		prefix = os.path.join(ffdir, 'R%dT%d_%s' % (rho, Temp,cgtype))
		log, hist = pp.parseFileNames(prefix)
		hist = pp.parseHist(hist, '', []) 
		hist = hist['LD']
		r_t = hist[0]; h_t = hist[1]; r_m = hist[2]; h_m = hist[3]
		ax2.plot(r_t+i, h_t, linestyle = 'solid', color = clrs[i])
		ax2.plot(r_m+i, h_m, linestyle = 'dashed', color = clrs[i])
		ax2.set_xlabel('Local Density')
		ax2.set_ylabel('LD Distribution')
		txt = '-   AA\n--  CG'
		ax2.annotate(txt, xy = (5,0.20))
		ax2.annotate('Pair spline\n+\nLD potential', xy = (18,0.20))
		ax1.yaxis.set_major_locator(tickloc)
		
		savedir = os.getcwd()
		figname = os.path.join(savedir, 'plot3.png') 
		plt.savefig(figname, bbox_inches = 'tight', pad_inches = 0., dpi = 300)


def plot4():
	''' plots per atom SASA for AA, SP and SPLD cases
	'''
	target_dir = os.path.abspath('./data/analysis/feb15_runs_fsw')
	trjtypes = ['AA_lj', 'CG_lj_SP', 'CG_lj_SPLD']
	
	setup()
	fig = plt.figure(figsize = (4,4))
	ax = fig.add_subplot(111)
	plt.subplots_adjust(left = 0.22, bottom = 0.22)
	
	clrs = ['red', 'blue', 'green']
	for i, trjtype in enumerate(trjtypes):
		datafile = os.path.join(target_dir, '%s_hist1D_SASA_atom.pickle' % trjtype)
		data = pickle.load(open(datafile, 'r'))
		x = data['bin_centers']; y = data['bin_vals']
		x = x[1:-2]; y = y[1:-2]
		lbl = trjtype.split('_')[0]
		if not lbl == 'AA':
			lbl += '_' + trjtype.split('_')[2]
		ax.plot(x,y, color = clrs[i], marker = 'o', markeredgecolor = 'black', markersize = 5,
				 label = lbl)
				 
	ax.set_xlabel('SASA/atom ' + r'$(\AA^2)$')
	ax.set_ylabel('SASA/atom  Distribution')
	ax.set_xticks([10, 14, 18, 22, 26, 30])
	ax.set_yticks([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
	ax.legend(loc = 'best', prop = {'size': 12})

	
def plot5():
	'''
	plots 2D pmf for Rg vs R_EE for AA, SP and SPLD cases	
	'''
	
	target_dir = os.path.abspath('./data/analysis/feb15_runs_fsw')
	trjtypes = ['AA_lj', 'CG_lj_SP', 'CG_lj_SPLD']
	
	setup()
	fig = plt.figure(figsize = (6.6,3))
	plt.subplots_adjust(bottom = 0.2, wspace = 0.001)
	
	ax1 = plt.subplot(1,3,1)
	ax2 = plt.subplot(1,3,2, sharey = ax1)
	ax3 = plt.subplot(1,3,3, sharey = ax2)
	axlist = [ax1, ax2, ax3]
	
	for i, trjtype in enumerate(trjtypes):
		ax = axlist[i]
		datafile = os.path.join(target_dir, '%s_hist2D_Rg_R_EE.pickle' % trjtype)
		data = pickle.load(open(datafile, 'r'))
		Rg = data['bin_centers'][0]
		R_EE  = data['bin_centers'][1]
		pmf = data['pmf']
		lbl = trjtype.split('_')[0]
		if not lbl == 'AA':
			lbl += ' ' + trjtype.split('_')[2]
		
		X,Y = np.meshgrid(Rg, R_EE)
		Z = pmf * conv_factor
		Z = np.clip(pmf, -1000000, 1)
		Z = (Z- Z.min()) / (Z.max() - Z.min())
		X = (X- X.min()) / (X.max() - X.min())
		Y = (Y- Y.min()) / (Y.max() - Y.min())	
		
		im = ax.imshow(Z, cmap = plt.get_cmap('jet'), extent = [X.min(),X.max(),Y.min(),Y.max()])
		if (i==2):
			pass
			#cbar = plt.colorbar(im, orientation = 'vertical', cmap = 'jet', 
			#		 			ticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
			#cbar.ax.set_xticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontsize = 11)
			
		ax.set_title(lbl, fontsize = 11)
		ax.set_xticks([0.2, 0.4, 0.6, 0.8])
		ax.set_yticks([0,0.2,0.4,0.6,0.8,1.])
		
		if (i==0):
			ax.set_ylabel(r'$(R_{EE} - R_{EE, min})/(R_{EE,max} - R_{EE,min})$', fontsize = 11)
		else:
			plt.setp(ax.get_yticklabels(), visible = False)
		if (i==1):
			ax.set_xlabel(r'$(Rg - Rg_{min})/(Rg_{max} - Rg_{min})$', fontsize = 11, weight = 'bold')
		

#plot1()
#plot2()
plot3()
#plot4()
#plot5()
plt.show()
