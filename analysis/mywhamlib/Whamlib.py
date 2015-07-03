#!/usr/bin/env python

import os, sys, pickle
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



class WHAM:
	
	######################################################################## 
	##04/14/15 - Unit conversion issues addressed:
	## Rules: All MD simulations use the Lammps 'real' system of units. 
	##		That means: energies are all in kcal, kcal/mol etc. Consequently
	##		beta passed to the specific WHAM workhorse routine should treat
	##		Boltzmann constant in kcal/mol-K. 
	##		Wham routines will always return non-dimensional pmf's 
	##		and that is how they will be stored. Any conversion should
	##		be done during post-processing such as plotting, comparing, etc
	########################################################################
							
	kB = 0.001982923700 # kcal/mol-K
	Temp = 298        # K
	beta = 1.0/(kB*Temp)
	invbeta = kB*Temp
	useParser = True
		
	
	def __init__(self, Prefix, TrajPrefix, trajdir, biasfile, Step, nbins):
		self.trajdir = trajdir
		data = np.loadtxt(biasfile)
		self.centers = data[:,0]
		self.FConst =  data[:,1]
		self.nbins = nbins
		self.Step  = Step
		self.Prefix = Prefix
		self.TrajPrefix = TrajPrefix	
		self.parse_PElog = None
		self.parse_Ulog = None
		self.histogram = None
		self.Method = None
		self.histfile = self.Prefix + '_wham_hist.pickle'	
		self.useParser = True
		self.doSubSample = False
		self.freq = 10
	
	def addParser(self, parsers = None):
		if not parsers:
			raise TypeError('Must supply parsers')
		else:
			print "Adding parsers"
			self.parse_PElog = parsers[0]
			self.parse_Ulog =  parsers[1]

	
	def subSample(self):
		# reduce number of frames to reduce any correlation
		print 'Subsampling data...'
		self.E_pe = self.E_pe[:, 0:self.E_pe.shape[1]:self.freq]
		self.E_u = self.E_u[:,0:self.E_u.shape[1]:self.freq]
		self.Rg = self.Rg[:, 0:self.Rg.shape[1]:self.freq] 	
		self.Step['nframes'] /= self.freq
		
			
	def genHist(self):	
		# check if data pickle already exists	
		if os.path.isfile(self.histfile):
			print 'Loading parsed data from whamhistpickle'
			f = open(self.histfile, 'r')
			(E_pe, E_u, Rg, centers, FConst, histogram) = pickle.load(f)
			self.E_pe = E_pe
			self.E_u = E_u
			self.Rg = Rg
			self.histogram = histogram
			self.Step['nframes'] = self.E_pe.shape[1] # reset nframes to prevent errors with subsampling
			return			
	
		# subsample data
		if self.doSubSample:
			self.subSample()		

		if self.useParser and (not self.parse_PElog	 or not self.parse_Ulog):
			raise ValueError('Must define parser.')
			exit()

		print "Initialising WHAM for " + self.Prefix + "..."
		f = open(self.histfile, 'w')
		
		# to use separately defined parsers they may not be taken from
		# from the Parsers lib in analysis/utils/Parsers.py but declared separately
		# using the format f(trajdir, trjprefix, Step = {nstages, nsteps})
		if self.useParser:
			self.E_pe = self.parse_PElog(self.trajdir, self.TrajPrefix, self.Step)
			(self.E_u, self.Rg) = self.parse_Ulog(self.trajdir, self.TrajPrefix, self.Step)
		

		# Setup hists with subsample data
		self.SetupHist(self.Rg)
		
		pickle.dump((self.E_pe, self.E_u, self.Rg, self.centers, self.FConst, self.histogram), f)
		f.close()
		

	
	def SetupHist(self, Rg):
		nframes = self.Step['nframes']
		nstages = self.Step["nstages"]
	
		count = []
		bin_assignment = np.zeros([nstages, nframes])
		Rg_min = Rg.min()
		Rg_max = Rg.max()
		bin_min = 0.98 * Rg_min; bin_max = 1.02*Rg_max
		bin_delta = (bin_max - bin_min)/float(self.nbins);
		bin_center = np.zeros([self.nbins], np.float64)
		
		for i in range(self.nbins):	
			bin_center[i] = Rg_min + bin_delta/2 + bin_delta*i
		
		for i in range(nstages):
			bin_val = np.zeros([self.nbins], np.float64)
			for n in range(nframes):		
				x = int((Rg[i,n] - bin_min)/bin_delta)
				bin_assignment[i,n] = x
				bin_val[bin_assignment[i,n]] += 1.00
			bin_val /= (np.sum(bin_val) * bin_delta)
			count.append(bin_val)
		count = np.array(count)
	
		hist = {'delta': bin_delta, 'bin_center': bin_center, 
			'bin_assignment': bin_assignment, 'count' : count,
			'nbins': self.nbins}
		
		self.histogram = hist
		
	
	def return_pmf(self, pmf):
		if self.Method == 'pymbar_wham':	
			x = pmf[0]
			y = pmf[1]
			z = pmf[2]
		if self.Method == 'grossfield_wham':
			x = pmf[:,0]
			y = pmf[:,1]
			z = []
		if self.Method == 'my_wham':
			x = pmf[0]
			y = pmf[1]
			z = []
		return (x,y,z)
	
	
	def compute_pmf(self):
		self.freefile = self.Prefix + "_" + self.Method + "_freefile.pickle"
		
		if os.path.isfile(self.freefile):
			f = open(self.freefile, 'r')
			pmf = pickle.load(f)
			return self.return_pmf(pmf)

		if not self.Method:
			raise ValueError('No WHAM method defined')
			exit()

		if not os.path.isfile(self.histfile):
			raise ValueError('Histograms missing. Cannot compute PMF')
			exit()
		
		Peng = {'pe': self.E_pe, 'u': self.E_u}
		Pumbrella  = {'Rg': self.Rg, 'Rg0' : self.centers, 'FConst' : self.FConst}
		
		print "Using WHAM routine: ", self.Method
		s = 'from %s import main as pmf_workhorse' % self.Method
		exec(s)
		pmf = pmf_workhorse(self.beta, Peng, Pumbrella, self.histogram, self.Step)	
		pickle.dump(pmf, open(self.freefile, 'w'))
		return self.return_pmf(pmf)
		
	
	def plot_overlap(self):		
		if not os.path.isfile(self.histfile):
			raise ValueError('Histograms missing')
			exit()
			
		count = self.histogram['count']
		bin_center = self.histogram['bin_center']
		delta = self.histogram['delta']
			
		fig = plt.figure()
		plt.xlabel(r'$Rg (\AA)$', fontsize = 20)
		plt.ylabel('Distribution', fontsize = 20)
		for i in range(self.Step['nstages']):
			plt.plot(bin_center, count[i], linestyle = "solid", linewidth = 2)
			for j in range(self.nbins):
				plt.plot(bin_center[j]*np.ones([100]), np.linspace(0,count[i][j],100), linestyle = "dashed", linewidth = 2, color = "black")
				plt.gca().add_patch(Rectangle(((bin_center[j] - delta/2), 0), delta, count[i][j], facecolor="grey", alpha = 0.1))	

		#plt.show()
		plt.savefig('%s_overlap.svg' % self.Prefix, bbox_inches = 'tight', dpi = 300)
		return fig


	def refineParams(self, pmf, prod):
		"""
		Uses trial and error based methods to give estimates of umbrella widths 
		and force-constant values for different umbrellas
		"""
		
		if not os.path.isfile(self.histfile):
			raise ValueError('Histograms missing')
			exit()
		

		## Unwrap data
		pmf = pmf[np.isfinite(pmf)]
		nstages = self.Step['nstages']
		nframes = self.Step['nframes']
		count = self.histogram['count']		
		bin_center = self.histogram['bin_center']
		delta0 = [self.centers[i+1] - self.centers[i] for i in range(len(self.centers)-1)]
		Rg0 = [self.centers[0]]
		k_u = [self.FConst[0]]
		
		prod = 1.0

		for i in range(nstages-1):
			this_sigma = np.std(count[i])
			next_sigma = np.std(count[i+1])
			slope = (pmf[i+1] - pmf[i]) / (bin_center[i+1] - bin_center[i])
			
			k_new = k_u[-1] * (this_sigma/next_sigma)**2
			delta_new = prod/(abs(slope) + 0.001)
			Rg0_new = np.arange(Rg0[-1],self.centers[i+1],delta_new)
	
			[Rg0.append(this) for this in Rg0_new[1:]]			
			[k_u.append(this) for this in k_new * np.ones([len(Rg0_new)])]
	
		self.refinefile = self.Prefix + "_refine.txt"
		dump = np.zeros([len(Rg0),2], np.float64)
		for i in range(len(Rg0)):
			dump[i,0] = Rg0[i]
			dump[i,1] = k_u[i]
		np.savetxt(self.refinefile, dump, header = 'Rg\t\t\t FConst')
