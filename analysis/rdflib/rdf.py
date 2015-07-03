#/usr/bin/env python

import os, sys, pickle
import numpy as np
import sim
from scipy.integrate import simps
import time

################## ATOMTYPE DEFINITIONS ################################
monomer_atomtype = 3
water_atomtype = 2 #oxygen-atom only
hydrogen_atomtype = 1
atomtypedefs = [monomer_atomtype, water_atomtype, hydrogen_atomtype]
########################################################################


#################### RECOMPILE BACKGROUND FORT-CODE ####################
libRecompile = True
if libRecompile:
	curr_dir = os.getcwd()
	os.chdir(os.path.expanduser('~/c25ld/analysis/rdflib'))
	s = '''
gfortran -c RDFlib.f90
f2py -c -m RDF_fortlib RDFlib.f90 --fcompiler=gfortran
'''
	os.system(s)
	os.chdir(curr_dir)
########################################################################



sys.path.append(os.path.expanduser('~/c25ld/analysis'))
import RDF_fortlib as fortlib

def getPickle(TrajName):
	'''
	Pickle the trajectory for faster access
	'''
	pickleName = TrajName + '.pickle'
	if os.path.isfile(pickleName):
		print 'Loading from pickle...'
		Trj = pickle.load(open(pickleName, 'r'))
	else:
		print 'Pickling trajectory...'
		Trj = sim.traj.lammps.Lammps(TrajName)
		pickle.dump(Trj, open(pickleName, 'w'))
	return Trj


## Not implemented in this code
def getMappedTraj(TrajName):
	''' Map the traj to oxygen and 
		monomer atoms only for **much** faster
		frame stepping
	'''
	
	MappedTrajName = TrajName.split('.lammpstrj.gz')[0] + '_mapped.lammpstrj.gz'
	if os.path.isfile(MappedTrajName):
		return getPickle(MappedTrajName)
		

	Trj = getPickle(TrajName)

	n_mon = 25
	n_water = 1700
	
	mappedatoms = []
	for (i, AtomType) in enumerate(Trj.AtomTypes):
		if AtomType == monomer_atomtype or AtomType == water_atomtype:
			mappedatoms.append(i)
	
	Map = sim.atommap.PosMap()
	for (i,j) in enumerate(mappedatoms):
		if i < n_water:
			mappedatomname =  water_atomtype
		else:
			mappedatomname = monomer_atomtype
		Map += [sim.atommap.AtomMap(Atoms1 = j, Atom2 = i, Atom2Name = mappedatomname)]
	AtomTypes = [water_atomtype] * len(mappedatoms) + [monomer_atomtype]*len(mappedatoms)
	MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes)  

	sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, MappedTrajName, Verbose = True)  
	
	del Trj
	Trj = getPickle(MappedTrajName)
	
	return Trj




def getRDF(TrajName, Nbins, Cut = 20.0, Prefix = 'rdf', target_dir = os.getcwd(), 
			frame_start = 0, frame_stop = 260000, stepfreq = 1):
	
	''' Makes a RDF from supplied sim trajectory
		object. Only considers atoms within specified
		Cut
	'''
	
	# Get un-mapped trajectory
	Trj = getPickle(TrajName)
	
	## Parse trajectory and make frame stepper
	init = Trj[0]
	BoxL = Trj.FrameData['BoxL']; BoxL = BoxL[0]
	#check for non-periodic box 05/01/2015
	if BoxL == 0:
		BoxL = float(raw_input('BoxL = 0, found. Enter nonperiodic boxlength: '))
	FrameRange = range(frame_start, frame_stop, stepfreq)
	NFrames = len(FrameRange)
	AtomTypes = Trj.AtomTypes
	
	## Make bins
	Bin_min = 1.00
	Bin_max = Cut
	Bin_delta = (Bin_max - Bin_min)/float(Nbins)
	Bin_centers = np.zeros([Nbins], np.float64)
	for i in range(Nbins):
		Bin_centers[i] = Bin_min + (i+0.5)*Bin_delta
		
	pb = sim.utility.ProgressBar(Text = 'Processing frame by frame...', Steps = NFrames)
	g = np.zeros([Nbins, 3], np.float64) #1-MM #2-WW #3-MW
	 
	## Frame stepping
	for frame in FrameRange:
		Pos = Trj[frame]
		
		# Call fortran routine to update histogram counts every frame
		g = fortlib.find_rdf(g = g, bin_centers = Bin_centers, bin_delta = Bin_delta, 
					   pos = Pos, atomtypes = AtomTypes, boxl = BoxL, atomtypedefs = atomtypedefs)
								
		
		pb.Update(int(frame/stepfreq))
		
	
	## Normalize the rdf and
	rdf = {'bin_centers': Bin_centers, 'MM': '', 'WW': '', 'MW': ''}
	N_mon = len(np.where(AtomTypes == monomer_atomtype)[0]); print N_mon
	N_water = len(np.where(AtomTypes == water_atomtype)[0]); print N_water
	
	print 'Normalizing g(r)...'
	BoxVol = BoxL**3.
	g /= NFrames
	gofr = g * BoxVol
	for j in range(Nbins):
		r = Bin_centers[j] - 0.5*Bin_delta
		next_r = Bin_centers[j] + 0.5*Bin_delta
		gofr[j,:] /= (4.*np.pi/3.) * (next_r**3. - r**3.)

	rdf['MM'] = {'g': gofr[:,0]/(N_mon * (N_mon - 1)/2.), 'count': g[:,0]}
	rdf['WW'] = {'g' :gofr[:,1]/(N_water * (N_water - 1)/2.), 'count': g[:,0]}
	rdf['MW'] = {'g': gofr[:,2]/(N_mon * N_water), 'count': g[:,0]}
		
		
	# pickling data
	pickleName = os.path.join(target_dir, Prefix + '.pickle')
	pickle.dump(rdf, open(pickleName, 'w'))	
	
	return rdf



def calcCoeff(Cut):
	'''
	Calculates Local Density Coefficients 
	based on given UpperCut. It is assumed that
	LowerCut = 0.8 * UpperCut
	'''
	
	R2sq = Cut*Cut
	R1sq = 0.8 * 0.8 * R2sq
	ratio = R1sq/R2sq
	denom = (1-ratio)*(1-ratio)*(1-ratio)
	
	c0 = (1.-3.*ratio)/denom
	c2 = (1./R2sq) * (6*ratio)/denom
	c4 = - (1./(R2sq*R2sq)) * (3. + 3.*ratio)/denom
	c6 = (1./(R2sq*R2sq*R2sq)) * (2./denom)
	
	return np.array([c0,c2,c4,c6])
	



def getLocalDensity(TrajName, LDAtomType, LDUpperCuts, Prefix = 'LD', target_dir = os.getcwd(), 
					frame_start = 0, frame_stop = 260000, stepfreq = 1):

	'''
	Calculates Local Density of central atomtype due to neighboring 
	atom type and for different Upper Cuts. Currently used for estimating
	Local Densities of coarse grained SPC/E water model
	'''
	
	# Get un-mapped trajectory
	Trj = getPickle(TrajName)
	

	## Parse trajectory and make frame stepper
	init = Trj[0]
	BoxL = Trj.FrameData['BoxL']; BoxL = BoxL[0]
	FrameRange = range(frame_start, frame_stop, stepfreq)
	NFrames = len(FrameRange)
	AtomTypes = Trj.AtomTypes
	
	## Make local density count array
	NCuts = len(LDUpperCuts)
	NAtom = len(AtomTypes)
	Central_atomtype = LDAtomType[0]
	Neigh_atomtype = LDAtomType[1]
	LD = np.zeros([NAtom, NCuts], np.float64)
	
	## Precompute LD coefficients
	coeff = np.zeros([NCuts, 4]) 
	for i, cut in enumerate(LDUpperCuts):
		coeff[i] = calcCoeff(cut)
		
	pb = sim.utility.ProgressBar(Text = 'Processing frame by frame...', Steps = NFrames)
	
	# Frame stepping
	for frame in FrameRange:
		Pos = Trj[frame]
		
		# Call fortran routine for per-frame calculation
		LD = fortlib.find_localdensity(ld = LD, pos = Pos, coeff = coeff, lduppercuts = LDUpperCuts, boxl = BoxL,
								 atomtypes = AtomTypes, central_atomtype = Central_atomtype, neigh_atomtype = Neigh_atomtype)
																			
		pb.Update(frame/stepfreq)
		
	## Average over all the frames
	LD /= NFrames
	
	# select only the oxygens
	LD_oxygen = np.zeros([NAtom/3, NCuts], np.float64)
	print 'Selecing oxygen Local Densities...'
	count = 0
	for i in range(0,NAtom-2,3):
		LD_oxygen[count, :] = LD[i,:]
		count += 1
	
	
	# pickling data
	pickleName = os.path.join(target_dir, Prefix + '.pickle')
	pickle.dump((LDUpperCuts, LD_oxygen), open(pickleName, 'w'))
	
	return LD





def getFirstShellWater(TrajName, LDUpperCuts, FirstShellCut = 6.0, Prefix = 'LD_wm', target_dir = os.getcwd(),
					   frame_start = 0, frame_stop = 100, stepfreq = 100):
	'''
	Calculates the local densities and first shell waters
	for each monomer for several different cutoffs. 
	Also generates the monomer local density distribution.
	'''
	# Get un-mapped trajectory
	Trj = getPickle(TrajName)
	

	## Parse trajectory and make frame stepper
	init = Trj[0]
	BoxL = Trj.FrameData['BoxL']; BoxL = BoxL[0]
	FrameRange = range(frame_start, frame_stop, stepfreq)
	NFrames = len(FrameRange)
	AtomTypes = Trj.AtomTypes
	
	## Make Local density count array
	NCuts = len(LDUpperCuts)
	N_mon = len(np.where(AtomTypes == 3)[0])
	N_water = len(np.where(AtomTypes == 2)[0])
	LD_mon = np.zeros([N_mon*NFrames,NCuts], np.float64)
	FirstShellWaters = np.zeros(N_mon*NFrames, np.float64)
	
	
	## Precompute LD coefficients
	coeff = np.zeros([NCuts, 4], np.float64)
	for i, cut in enumerate(LDUpperCuts):
		coeff[i] = calcCoeff(cut)
			
	pb = sim.utility.ProgressBar(Text = 'Processing frame by frame...', Steps = NFrames)

	# Frame stepping
	count = 0
	for frame in FrameRange:
		Pos = Trj[frame]
		
		# Call fortran subroutine for per-frame computation
		fsw = np.zeros(N_mon, np.float64)
		ld = np.zeros([N_mon, NCuts], np.float64)
		(ld, fsw) = fortlib.find_firstshellwaters(ld = ld, fsw = fsw, pos = Pos, atomtypes = AtomTypes, n_water = N_water,
												 lduppercuts = LDUpperCuts, firstshellcut = FirstShellCut, 
												 coeff = coeff, boxl = BoxL, atomtypedefs = atomtypedefs)
		
		## Average 
		#LD_mon += ld
		#FirstShellWaters += fsw										 										 
		
		## Per-frame
		FirstShellWaters[count:count+N_mon] = fsw
		for j in range(NCuts):
			LD_mon[count:count+N_mon,j] = ld[:,j]
		count += N_mon
		
		pb.Update(frame/stepfreq)
		
	
	#LD_mon /= NFrames
	#FirstShellWaters /= NFrames
	
	# pickling data
	pickleName = os.path.join(target_dir, Prefix + '.pickle')
	pickle.dump((LDUpperCuts, LD_mon, FirstShellWaters), open(pickleName, 'w'))
		
	return (FirstShellWaters, LD_mon)	
