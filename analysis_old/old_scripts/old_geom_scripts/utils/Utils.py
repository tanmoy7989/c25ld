#/usr/bin/env python

import sim
import numpy as np
import os, sys, pickle

from sim.utility import ProgressBar as pb

sys.path.append(os.path.expanduser('~/c25ld'))
from chem_data import getModelData
d = getModelData()
poly_dict = d['poly_dict']

MinImageDEBUG = False


class Geom:

		''' Class for defining computation for geometrical properties
		of trajectories.Defines compute routines for all-frame based 
		trajectory paramters. 
		'''
		#TODO: Write the looping in fortran
		# 04/14/15 : Now, can handle multiple polymer chains
	
		isDataPickled = False	
		isPickled = False
		doMinimage = True
		
		hasRg = True
		hasR_EE = True
		hasKappa = True
		hasSASA = True
		N_poly = 1	#default
		N_mon = 25	#default
		

		def __init__(self, traj, Prefix = 'geom', target_dir = os.getcwd(), useBoxL = ''):
			self.traj = traj
			self.data = os.path.join(os.path.abspath(target_dir), Prefix + '.pickle')
			self.useBoxL = useBoxL
			
			if os.path.isfile(self.data):
				print self.data
				self.isDataPickled = True
				self.isPickled = True
			else:
				if os.path.isfile(self.traj + '.pickle'):
					self.isPickled = True
				else: 
					self.pickleTrj()

		def pickleTrj(self):
			pickleName = self.traj + '.pickle'
			print "Pickling trajectory..."
			Trj = sim.traj.lammps.Lammps(self.traj)
			pickle.dump(Trj, open(pickleName, 'w'))
			self.isPickled = True


		def compute(self, groupID):
			
			if self.isDataPickled:
				print "Data pickle already present. Loading ..."
				data = pickle.load(open(self.data, 'r'))
				return data

			if not self.isPickled:
				raise TypeError('Reading unpickled trajectory is too time consuming')


			pickleName = self.traj + '.pickle'
			Trj = pickle.load(open(pickleName, 'r'))
			init = Trj[0]
			if not self.useBoxL:
				BoxL = Trj.FrameData['BoxL']
			else:
				BoxL = np.array([1.,1.,1.]) * self.useBoxL
			
			ind = np.where(Trj.AtomTypes == groupID)[0]
			Nframes = len(Trj)
			computes = {}
			if self.hasRg:
				Rg = np.zeros(Nframes)
				computes['Rg'] = np.zeros(Nframes)
			if self.hasR_EE:
				R_EE = np.zeros(Nframes)
				computes['R_EE'] = np.zeros(Nframes)
			if self.hasKappa:
				l = np.zeros([Nframes,3])
				kappa = np.zeros(Nframes)
				computes['lambda'] = np.zeros([Nframes, 3])
				computes['kappa'] = np.zeros(Nframes)
			if self.hasSASA:
				Radius = 0.5 * poly_dict['LJSigma'] * 2**(1./6.)
				Radius = 0.5*(1.4 + Radius)
				SASA_atom = np.zeros([Nframes, self.N_mon])
				SASA = np.zeros(Nframes)
				computes['SASA_atom'] = np.zeros([Nframes, self.N_mon])
				computes['SASA'] = np.zeros(Nframes)

			frameStatus = pb(Text = 'Reading frame by frame', Steps = Nframes)
			for frame in range(Nframes):
				Pos = Trj[frame][ind]
				# Reimaging current frame
				if self.doMinimage:
					for i in range(1, self.N_mon * self.N_poly):
    						Pos[i] = sim.geom.Reimage(Pos[i], Pos[i-1], BoxL)
				if MinImageDEBUG:
				    print BoxL
				    print Pos; raw_input()
				## per-polymer
				for j in range(self.N_poly):
					Pos_poly = Pos[j*self.N_mon : (j+1)*self.N_mon]
					
					## per-frame based geometrical computations
					if self.hasKappa:
						this_kappa = 0.
						this_l = np.zeros(3, np.float64)
						diff = Pos_poly - Pos_poly.mean(axis = 0)
						for ax in range(3):
							this_l[ax] = np.sqrt(np.sum(diff[:,ax]**2, axis = 0))
						this_kappa = (1.5) * (this_l[0]**4 + this_l[1]**4 + this_l[2]**4) / (this_l[0]**2 + this_l[1]**2 + this_l[2]**2)**2
						this_kappa -= 0.5	
						l[frame, :] += this_l
						kappa[frame] += this_kappa
					
					if self.hasRg:
						this_Rg = 0
						this_Rg = np.sqrt(np.sum((Pos_poly-Pos_poly.mean(axis=0))**2)/self.N_mon) 
						Rg[frame] += this_Rg
						
					if self.hasR_EE:
						this_R_EE = 0.
						this_R_EE = np.sqrt(np.sum((Pos_poly[-1]-Pos_poly[0])**2))
						R_EE[frame] += this_R_EE	
						
					if self.hasSASA:
						this_SASA_atom = np.zeros(self.N_mon, np.float64)
						this_SASA = 0.
						Radii = Radius * np.ones([self.N_mon])
						this_SASA_atom = sim.geom.SphereSurfaceAreas(Pos_poly, Radii)
						this_SASA = np.sum(this_SASA_atom)
						SASA_atom[frame,:] += this_SASA_atom 
						SASA[frame] += this_SASA

					
				l[frame,:] /= self.N_poly
				kappa[frame] /= self.N_poly
				SASA_atom[frame, :] /= self.N_poly
				SASA[frame] /= self.N_poly
				Rg[frame] /= self.N_poly
				R_EE[frame] /= self.N_poly
				
				frameStatus.Update(frame)
				
			if self.hasRg:
				computes['Rg'] = Rg
			if self.hasR_EE:
				computes['R_EE'] += R_EE
			if self.hasKappa:
				computes['lambda'] += l
				computes['kappa'] += kappa
			if self.hasSASA:
				computes['SASA_atom'] += SASA_atom
				computes['SASA'] += SASA
	
			pickle.dump(computes, open(self.data, 'w'))
			self.isDataPickled = True

			return computes


class Hist:

	'''
	Class for generating 1D and 2D histograms.
	Can optionally also calculate PMF by inverting said histograms
	'''

	isDataPickled = False

	def __init__(self, z, nbins, Prefix, target_dir = os.getcwd(), Dim = 1):
		
		self.data = os.path.join(os.path.abspath(target_dir), Prefix + '.pickle')
		if os.path.isfile(self.data):
			self.isDataPickled = True

		self.nbins = nbins
		self.z = z
		self.Dim = Dim
		if not self.Dim == 1:
			if not len(self.z)>1 or not len(self.nbins)>1:
				raise TypeError('Cannot do 2D histograms with 1 variable')


	def makeHist(self, normalize = True, doPMF = True):
		if self.isDataPickled:
			return

		if not self.Dim == 1:
			raise TypeError('Variable # mismatch')

		z = self.z
		Nframes = len(z)
		bin_min = 0.98 * z.min(); bin_max = 1.02*z.max()
		delta = (bin_max - bin_min)/float(self.nbins)
		bin_centers = np.zeros(self.nbins)
		bin_vals = np.zeros(self.nbins)
		pmf = np.zeros(self.nbins)
		for i in range(self.nbins):
			bin_centers[i] = bin_min + (i+0.5) * delta
			
		frameStatus = pb(Text = 'Binning frame by frame', Steps = Nframes)
		for i in range(Nframes):
		
			assignment = int((z[i] - bin_min)/delta)
			bin_vals[assignment] += 1.0
		
			frameStatus.Update(i)
		
		if normalize:
			#bin_vals /= (np.sum(bin_vals) * delta)
			bin_vals /= np.trapz(bin_vals, bin_centers, dx = delta)
		if doPMF:
			pmf = - np.log(bin_vals)
		

		hist = {'bin_centers': bin_centers, 'bin_vals': bin_vals, 'pmf' : pmf}
		pickle.dump(hist, open(self.data, 'w'))
		self.isDataPickled = True


	def makeHist2D(self, normalize = True, doPMF = True):
		if self.isDataPickled:
			return

		if not self.Dim == 2:
			raise TypeError('Variable # mismatch')

		nbins1 = self.nbins[0]; nbins2 = self.nbins[1]	
		z1 = self.z[0]; z2 = self.z[1]
		if not len(z1)==len(z2):
			raise TypeError('Unequal counts of z1 and z2')
		
		Nframes = len(z1)
		bin1_min = 0.98*z1.min(); bin1_max = 1.02*z1.max()
		bin2_min = 0.98*z2.min(); bin2_max = 1.02*z2.max()
		delta1 = (bin1_max - bin1_min)/float(nbins1); delta2 = (bin2_max - bin2_min)/float(nbins2)
		
		bin_centers_1 = np.zeros(nbins1); bin_centers_2 = np.zeros(nbins2)
		bin_vals = np.zeros([nbins1, nbins2])
		pmf = np.zeros([nbins1, nbins2])
		
		for i in range(nbins1):
			bin_centers_1[i] = z1.min() + (i+0.5) * delta1
		for i in range(nbins2):
			bin_centers_2[i] = z2.min() + (i+0.5) * delta2

		frameStatus = pb(Text = 'Reading frame by frame', Steps = Nframes)
		for i in range(Nframes):
				
			assignment1 = int((z1[i] - bin1_min)/delta1)
			assignment2 = int((z2[i] - bin2_min)/delta2)

			bin_vals[assignment1, assignment2] += 1.0
			
			frameStatus.Update(i)

		if normalize:
			bin_vals /= (np.sum(bin_vals) * delta1 * delta2)

		if doPMF:
			pmf = - np.log(bin_vals)

		hist = {'bin_centers': [bin_centers_1, bin_centers_2], 'bin_vals': bin_vals, 'pmf' : pmf}
		pickle.dump(hist, open(self.data, 'w'))
		self.isDataPickled = True


