import os, sys, pickle
import numpy as np
import sim

sys.path.append(os.path.abspath('./rdflib'))
import RDF_fortlib as fortlib

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


## Parse user input
TrajName = sys.argv[1]
target_dir = sys.argv[2]
Prefix = sys.argv[3]

## Other parameters
#if not Prefix:
#	Prefix = 'firstshell_spce'
FirstShellCut = 3.0
LDUpperCuts = [3.0, 3.15, 3.18, 3.2, 3.22, 3.25, 3.35, 3.5, 3.6, 3.75, 4.0, 5.0, 6., 8.0]
frame_start = 0
frame_stop = 20000
stepfreq = 100

# Get un-mapped trajectory
pickleName = TrajName + '.pickle'
if not os.path.isfile(pickleName):
	print 'Pickling data...'
	Trj = sim.traj.lammps.Lammps(TrajName)
	pickle.dump(Trj, open(pickleName, 'r'))
else:
	Trj = pickle.load(open(pickleName, 'r'))
	
## Parse trajectory and make frame stepper
init = Trj[0]
BoxL = Trj.FrameData['BoxL']; BoxL = BoxL[0]
FrameRange = range(frame_start, frame_stop, stepfreq)
NFrames = len(FrameRange)
AtomTypes = Trj.AtomTypes
	
## Make Local density count array
NCuts = len(LDUpperCuts)
NAtom = len(AtomTypes)
NWater = NAtom/3
LD = np.zeros([NWater*NFrames,NCuts], np.float64)
FSW = np.zeros(NWater*NFrames, np.float64)
		
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
	fsw = np.zeros(NAtom, np.float64)
	ld = np.zeros([NAtom, NCuts], np.float64)
	(ld, fsw) = fortlib.find_spce_fsw(ld = ld, fsw = fsw, pos = Pos, atomtypes = AtomTypes,
												 lduppercuts = LDUpperCuts, firstshellcut = FirstShellCut, 
												 coeff = coeff, boxl = BoxL)
		
	# pack per-frame data into global arrays
	count_o2 = 0
	fsw_o2 = np.zeros(NWater, np.float64)
	ld_o2 = np.zeros([NWater,NCuts], np.float64)
	for i in range(0, NAtom,3):
		fsw_o2[count_o2] = fsw[i]
		ld_o2[count_o2,:] = ld[i,:]
		count_o2 += 1
	
	FSW[count:count+NWater] = fsw_o2
	for j in range(NCuts):
		LD[count:count+NWater,j] = ld_o2[:,j]
	count += NWater
		
	pb.Update(frame/stepfreq)
		
# pickling data
pickleName = os.path.join(target_dir, Prefix + '.pickle')
pickle.dump((LDUpperCuts, LD, FSW), open(pickleName, 'w'))			
