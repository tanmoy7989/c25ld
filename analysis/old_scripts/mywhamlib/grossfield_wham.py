import numpy as np
import os
import sys
sys.dont_write_bytecode = True

###############################################################################
## Issued with units -- fixed 04/16/15
## The file wham.h maintains the global unit system
## used for all the calculations. Since all Lammps simulation
## data are in the 'real' system of units (i.e. energies in kcal/mol)
## wham.h maintains the Boltzmann constant in kcal/mol-K
##############################################################################

## Unit conversions for returning non-dimensional pmf
kB = 0.001982923700 # kcal/mol-K
Temp = 298        # K
beta = 1./ (kB * Temp)

whamexec = os.path.expanduser('~/c25ld/analysis/mywhamlib/grossfield_wham/wham/./wham')
doErrAnalysis = True

def main(beta, eng, umbrella, hist, step):
	# Unwrapping dicts passed to the function
	nstages = step['nstages']
	nframes = step['nframes']
	
	E_pe = eng['pe']
	E_u  = eng['u']
	nbins = hist['nbins']
	Rg0 = umbrella['Rg0']
	FConst = umbrella['FConst']
	# Grossfield's routine works with V(x) = 1/2 k (x-x0)**2
	# similar to the user-colvars package in lammps
	
	Rg = umbrella['Rg']
	# Setting parameters
	Temp = 298.0
	hist_min = hist['bin_center'].min() - hist['delta']/2.
	hist_max = hist['bin_center'].max() + hist['delta']/2.
	TOL = 1e-5
	numpad = 0
	metadatafile = 'wham.mdata'
	freefile = "PMF.txt"
	
	n_MCtrials = 6
	seed = 82717
	
	# Writing meta-data file
	WriteMetadata(Rg, Rg0, E_u, FConst, nstages, nframes, Temp)
	cmdstring = "%g %g %d %g %g %d %s %s" % (hist_min, hist_max, nbins, TOL, Temp, numpad, metadatafile, freefile)
	if doErrAnalysis:
		cmdstring += " %d %d" % (n_MCtrials, seed)
	
	# Running Grossfield wham
	os.system(whamexec + " " + cmdstring)
	print 'Wham run complete...'
	pmf = np.loadtxt(freefile)
	Cleanup(nstages)

	# Return non-dimensional pmf
	pmf *= beta
	return pmf 
	
	
def WriteMetadata(Rg, Rg0, E_u, FConst, nstages, nframes, Temp):

	mdatastring = ''	
	
	# writing timeseries files
	data = np.zeros([nframes, 3], np.float64)
	print "Writing timeseries files..."
	for stage in range(nstages):
		trgfile = "tseriesfile" + str(stage) + ".txt"	
		for frame in range(nframes):
			data[frame, 0] = frame
			data[frame, 1] = Rg[stage,frame]
			data[frame, 2] = E_u[stage, frame] 
		
		np.savetxt(trgfile, data)
		## header addition supported only for later versions of python
		#np.savetxt(trgfile, data, header = "Timeseries file for stage " + str(stage))

		#index = np.where(E_u[stage,:] == np.min(E_u[stage,:]))[0][0]
		loc_min_bias = Rg0[stage]
		spring = FConst[stage]
		mdatastring += "%s %f %f" % (trgfile, loc_min_bias, spring)
		if doErrAnalysis:
			mdatastring += " %d %f" % (1, Temp)
		mdatastring += "\n"		
	
	# writing metadatafile
	print "Writing Meta-data-file"
	with open("wham.mdata", "w") as of:
		of.write("#Meta-data file for WHAM runs using Grossfield WHAM code on c25T298 system\n\n")
		of.write(mdatastring)
	of.close()



def Cleanup(nstages):
	for stage in range(nstages):
		seriesfile = "tseriesfile" + str(stage) + ".txt"
		if os.path.isfile(seriesfile):
			os.remove(seriesfile)
	if os.path.isfile("wham.mdata"): 
		os.remove("wham.mdata")

	if os.path.isfile("PMF.txt"):
		os.remove("PMF.txt")

