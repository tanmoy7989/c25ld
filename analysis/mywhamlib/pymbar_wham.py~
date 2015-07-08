import numpy as np
import pymbar
import timeseries

import sys
sys.dont_write_bytecode = True

################################################################################
## 04/14/15: Unit conversion issues addressed 
## Units passed from calling parent functions are as:
## Energy: kcal/mol, Rg: Angstroms, Force-const: kcal/mol-Angstrom^2, Temp: Kelvin
## pmf returned is dimensionless (in units of kT)
################################################################################


def main(beta, eng, umbrella, hist, step):

	# Unwrapping dicts passed to the function
	nstages = step['nstages']
	nframes = step['nframes']
	
	E_pe = eng['pe']
	E_u  = eng['u']
	u_kn = E_pe - E_u
	Rg_kn = umbrella['Rg']
	Rg0 = umbrella['Rg0']
	FConst = umbrella['FConst']
	
	bin_kn = hist['bin_assignment']
	nbins  = hist['nbins']
	bin_center = hist['bin_center']
	count = hist['count']
	N_k = np.ones(nstages, np.int32)*nframes
			
	# set minimum of u_kn to 0 : Arbitrary
	u_kn -= u_kn.min()
											
	# Evaluate correlated reduced energies in all umbrellas
	print "Evaluating reduced potential energies..."
	N_max = np.max(N_k)
	u_kln = np.zeros([nstages,nstages,N_max], np.float64)
	for k in range(nstages):
		for n in range(N_k[k]):
			u_kln[k,:,n] = (u_kn[k,n] + E_u[:,n]) * beta


	# Initialising MBAR
	print "Building MBAR object. May take some time..."
	mbar = pymbar.MBAR.MBAR(u_kln, N_k, verbose = True, method = 'Newton-Raphson')


	# Extracting PMF (in units of kT)
	print "Computing non-dimensional PMF..."
	#(PMF, dPMF) = mbar.computePMF(u_kn, bin_kn, nbins)
	
	## The above function is broken. So, below is my own version of calcualating 
	## the pmf. It assumes that N_k[k] = nframes for all k 
	## 04/16/15
	A = mbar.f_k
	PMF = np.zeros([nbins], np.float64)
	for i in range(nbins):
		U_bias_bin = 0.5 * FConst * (bin_center[i] - Rg0)**2.
		U_bias_bin *= beta
		PMF[i] = -np.log(np.sum(count[:,i])) + np.log(nframes) + logSumExp((A - U_bias_bin))
	
	return (bin_center, PMF, A)


def logSumExp(w):
	w_max = max(w)
	return w_max + np.log(np.sum(np.exp(w - w_max)))
	
	
