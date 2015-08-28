#/usr/bin/env python

import os, sys
import numpy as np
import pickle

sys.path.append(os.path.expanduser('~/c25ld/analysis'))
from utils import Utils
import sim


class Analysis:
	hist1D_flags = {'kappa' : True, 'R_EE': True, 'SASA' : True, 'Rg': True}

	def __init__(self, trjname, trjtype = 'AA', trjsrc = 'lmp', target_dir = os.getcwd(), n_mon = 25):

		self.n_mon = n_mon
		self.trj = trjname; 
		self.type = trjtype;
		self.src_code = trjsrc
		self.target_dir = os.path.abspath(target_dir)
		self.gPrefix = self.type + '_geom'
		self.h1D_Prefix = {}
		for key in self.hist1D_flags.keys():
			self.h1D_Prefix[key] = self.type + '_hist1D_' + key

		if self.src_code == 'lmp':
			self.groupID = 3
		else:
			self.groupID = 0

		self.data = {}

	def compute_geom(self):
		
		proplist = []
		[proplist.append(key) for key in self.hist1D_flags.keys() if self.hist1D_flags[key]]
		g = Utils.Geom(traj = self.trj, Prefix = self.gPrefix, target_dir = self.target_dir)
		## generalized chain length 
		## added by T.S. on 08/28/2015
		g.N_poly = 1
		g.N_mon = self.n_mon
		for key in self.hist1D_flags.keys():
			attr = 'has' + key
			if hasattr(g,attr):
				setattr(g,attr,self.hist1D_flags[key])
		
		print 'Computing %s for %s\n' % (str(proplist), self.gPrefix)
		self.data = g.compute(groupID = self.groupID) 

	def make_Hist1D(self,nbins):
		if not self.data:
			raise ValueError('Need to run compute_geom() method before setting up histograms')
		for key in self.hist1D_flags.keys():
			if self.hist1D_flags[key]:
				print 'Setting up histograms for ', self.h1D_Prefix[key]
				h = Utils.Hist(z = self.data[key], nbins = nbins, Prefix = self.h1D_Prefix[key], target_dir = self.target_dir, Dim = 1)
				h.makeHist(normalize = True)


	
	def make_SASA_atom(self, nbins):
		if not self.data:
			raise ValueError('Need to run compute_geom() method before setting up histograms')

		hPrefix = self.type + '_hist1D_' + 'SASA_atom'
		print 'Setting up histograms for ', hPrefix
		histPickle = os.path.join(self.target_dir, hPrefix + '.pickle')
		if os.path.isfile(histPickle):
			return
		
		Nframes = len(self.data['SASA_atom'])
		bin_centers = np.zeros(nbins)
		bin_vals = np.zeros(nbins)
		all_bin_vals = []
		pmf = np.zeros(nbins)
		
		pb = sim.utility.ProgressBar(Text = 'Binning frame..', Steps = Nframes)
		for frame in range(Nframes):
			tmp_Prefix = hPrefix + str(frame)
			z = self.data['SASA_atom'][frame]
			h = Utils.Hist(z = z, nbins = nbins, Prefix = tmp_Prefix, target_dir = self.target_dir, Dim = 1)
			h.makeHist(normalize = True)

			histPickle = os.path.join(self.target_dir, tmp_Prefix+'.pickle')
			this_hist = pickle.load(open(histPickle, 'r'))
			if frame == 0:
				bin_centers = this_hist['bin_centers']
			bin_vals += this_hist['bin_vals']
			all_bin_vals.append(bin_vals)
			os.remove(histPickle)
			pb.Update(frame)

		bin_vals /= Nframes
		pmf = - np.log(bin_vals)
		hist = {'bin_centers': bin_centers, 'bin_vals': bin_vals, 'pmf': pmf, 'all_bin_vals': all_bin_vals}
		histPickle = os.path.join(self.target_dir, hPrefix + '.pickle')
		pickle.dump(hist, open(histPickle,'w'))



	#TODO: Convert Rg into Rg binned from umbrella sampling
	def make_Rg_R_EE(self,nbins):
		if not self.data:
			raise ValueError('Need to run compute_geom() method before setting up histograms')
		hPrefix = self.type + '_hist2D_' + 'Rg_R_EE'
		print 'Setting up 2D histogram from Rg vs R_EE for ',hPrefix
		histPickle = os.path.join(self.target_dir, hPrefix + '.pickle')
		if os.path.isfile(histPickle):
			return
		h = Utils.Hist(z = [self.data['Rg'], self.data['R_EE']], nbins = nbins, Prefix = hPrefix, target_dir = self.target_dir, Dim = 2)
		h.makeHist2D(normalize = True)
	

##### MAIN #####

inputfile = sys.argv[1]

# Parse input
argdict = eval(file(inputfile).read())
trjname = argdict['name']
trjtype = argdict['type']
trjsrc = argdict['src']
target_dir = argdict['target_dir']
n_mon = argdict['n_mon']

a = Analysis(trjname = trjname, trjtype = trjtype, trjsrc = trjsrc, target_dir = target_dir, n_mon = n_mon)

a.compute_geom()
a.make_Hist1D(nbins = 50)
a.make_SASA_atom(nbins = 50)
a.make_Rg_R_EE(nbins = [100,100])
